from align.primitive.default.canvas import DefaultCanvas
from align.cell_fabric.generators import Region, Wire, Via
from align.cell_fabric.grid import EnclosureGrid, UncoloredCenterLineGrid, SingleGrid, CenteredGrid, CenterLineGrid
from math import floor, ceil
import collections
import string
import logging
logger = logging.getLogger(__name__)


class MOSGenerator(DefaultCanvas):

    def __init__(self, pdk, height, fin, gate, gateDummy, shared_diff, stack, bodyswitch, **kwargs):
        self.primitive_constraints = kwargs.get('primitive_constraints', [])
        self.primitive_parameters = kwargs.get('primitive_parameters')
        super().__init__(pdk)

        exact_width = None
        exact_length = None
        length_diff = 0
        dynamic_space = 0
        if self.primitive_parameters:
            device_names = [*self.primitive_parameters.keys()]
            exact_width = self.primitive_parameters[device_names[0]]['W']
            exact_width = int(float(exact_width)*1E9)  ### Width in nanometers
            exact_length = self.primitive_parameters[device_names[0]]['L']
            exact_length = round(float(exact_length)*1E9)  ### Length in nanometers
            assert exact_length % 2 ==0, f"Transistor gate length {exact_length} must be even"
            assert exact_width % 2 ==0, f"Transistor width {exact_width} must be even"
            length_diff = exact_length - self.pdk['Poly']['Width']
            assert length_diff >= 0, f"Transistor gate length {exact_length} must be greater than the minimum gate length {self.pdk['Poly']['Width']}"
            dynamic_space = 20 if length_diff >= 100 else 0 ### This changes poly pitch based on Lg dependent DRCs
            self.pdk['Poly']['Width'] = exact_length
            self.pdk['Poly']['Pitch'] = self.pdk['Poly']['Pitch'] + length_diff + dynamic_space
            self.pdk['Poly']['Offset'] = self.pdk['Poly']['Pitch']//2
            self.pdk['M1']['Pitch'] = self.pdk['Poly']['Pitch']

        self.exact_patterns = None
        # For planar CMOS: fin parameter represents width in number of fin pitches
        # Use exact_width if provided, otherwise calculate from fin parameter
        if exact_width:
            # Calculate equivalent fin count from exact width
            fin_equiv = max(4, 4*ceil((exact_width/self.pdk['Fin']['Pitch'])/4))
        else:
            fin_equiv = fin if fin else 4

        for const in self.primitive_constraints:
            if const.constraint == 'Generator':
                if const.parameters is not None:
                    if const.parameters.get('exact_patterns'):
                        self.exact_patterns = const.parameters.get('exact_patterns')
                        self.exact_patterns = self.exact_patterns[0]
                    if const.parameters.get('height'):
                        height = const.parameters.get('height')

        height = max(height, 4*ceil((fin_equiv+16)/4))
        self.finsPerUnitCell = height
        assert self.finsPerUnitCell % 4 == 0
        assert (self.finsPerUnitCell*self.pdk['Fin']['Pitch'])%self.pdk['M2']['Pitch']==0
        self.m2PerUnitCell = (self.finsPerUnitCell*self.pdk['Fin']['Pitch'])//self.pdk['M2']['Pitch']
        self.unitCellHeight = self.m2PerUnitCell* self.pdk['M2']['Pitch']
        ######### Derived Parameters ############
        self.shared_diff = shared_diff
        self.stack = stack
        self.bodyswitch = bodyswitch
        self.gateDummy = gateDummy
        self.gate = (2*gate)*self.stack
        self.gatesPerUnitCell = self.gate + 2*self.gateDummy*(1-self.shared_diff)
        # For planar CMOS, use exact_width to determine actual fin usage
        if exact_width:
            self.finDummy = (self.finsPerUnitCell - fin_equiv)//2
        else:
            self.finDummy = (self.finsPerUnitCell - fin)//2 if fin else (self.finsPerUnitCell - 4)//2
        self.lFin = height ### This defines numebr of fins for tap cells; Should we define it in the layers.json?
        assert self.finDummy >= 8, f"number of fins/width {fin_equiv if exact_width else fin} in the transistor must be less than unit cell height {self.finsPerUnitCell -2*8}"
        # Removed fin > 1 assertion for planar CMOS compatibility
        assert gateDummy > 0
        unitCellLength = self.gatesPerUnitCell* self.pdk['Poly']['Pitch']
        activeOffset = self.unitCellHeight//2 -self.pdk['Fin']['Pitch']//2
        activeWidth =  exact_width if exact_width else self.pdk['Fin']['Pitch']*fin
        activePitch = self.unitCellHeight
        RVTWidth = activeWidth + 2*self.pdk['Active']['active_enclosure']

        # Calculate Nselect width to properly enclose Active
        try:
            NselectEnclosure = self.pdk['Nselect'].get('Nselect_enclosure_Active', 60)
        except (KeyError, TypeError):
            NselectEnclosure = 60
        NselectWidth = activeWidth + 2 * NselectEnclosure

        # Store for use in addNMOSArray to calculate optimized NP region
        self.activeOffset = activeOffset
        self.activeWidth = activeWidth
        self.NselectEnclosure = NselectEnclosure

        # Optimize Poly extension:
        # - Bottom: minimum DRC requirement (60nm beyond Active)
        # - Top: align with Pc layer bottom
        poly_extension = self.pdk['Active'].get('PolyExtension', 60)

        # Calculate target positions
        active_bottom = activeOffset - activeWidth//2
        poly_bottom_target = active_bottom - poly_extension

        # Calculate Pc bottom position
        # Pc is at grid (m2PerUnitCell - 4) with offset = M2.Pitch
        pc_grid_index = self.m2PerUnitCell - 4
        pc_center = self.pdk['M2']['Pitch'] + pc_grid_index * self.pdk['M2']['Pitch']
        pc_width = self.pdk['Pc'].get('PcWidth', 100)
        poly_top_target = pc_center - pc_width // 2  # Pc bottom

        # Calculate offset and stoppoint for EnclosureGrid
        # poly_bottom = offset + stoppoint
        # poly_top = offset + pitch - stoppoint
        pitch = self.unitCellHeight
        offset = (poly_bottom_target + poly_top_target - pitch) // 2
        stoppoint = poly_bottom_target - offset

        self.pl = self.addGen( Wire( 'pl', 'Poly', 'v',
                                     clg=UncoloredCenterLineGrid( pitch= self.pdk['Poly']['Pitch'], width= self.pdk['Poly']['Width'], offset= self.pdk['Poly']['Offset']),
                                     spg=EnclosureGrid( pitch=pitch, offset=offset, stoppoint=stoppoint, check=False)))

        # Calculate stoppoint and align to M2.Pitch to ensure cell height is grid-aligned
        raw_stoppoint = self.pdk['V1']['VencA_L'] + self.pdk['M2']['Width']//2
        aligned_stoppoint = raw_stoppoint  # Use raw value to minimize M1 extension

        # Calculate M1 width to ensure proper V1 enclosure
        # M1 width must be >= V1 width + 2 * V1 enclosure by M1 requirement
        v1_width = self.pdk['V1']['WidthX']  # 100nm
        v1_enc_m1 = self.pdk['V1']['VencA_L']  # 40nm
        m1_width_for_v1 = v1_width + 2 * v1_enc_m1  # 180nm
        m1_width = max(self.pdk['M1']['Width'], m1_width_for_v1)

        self.m1_updated = self.addGen( Wire( 'm1_updated ', 'M1', 'v',
                                     clg=UncoloredCenterLineGrid( pitch=self.pdk['M1']['Pitch'], width=m1_width),
                                     spg=EnclosureGrid( pitch=self.pdk['M2']['Pitch'], stoppoint=aligned_stoppoint, check=False)))


        # Calculate M2 width to ensure proper V1 enclosure
        # M2 width must be >= V1 width + 2 * V1 enclosure by M2 requirement
        v1_height = self.pdk['V1']['WidthY']  # 100nm
        v1_enc_m2 = self.pdk['V1']['VencP_L']  # 40nm
        m2_width_for_v1 = v1_height + 2 * v1_enc_m2  # 180nm
        m2_width = max(self.pdk['M2']['Width'], m2_width_for_v1)
        self.m2_updated = self.addGen( Wire( 'm2_updated ', 'M2', 'h',
                                     clg=UncoloredCenterLineGrid( pitch=self.pdk['M2']['Pitch'], width=m2_width),
                                     spg=EnclosureGrid( pitch=self.pdk['M1']['Pitch'], stoppoint=self.pdk['V1']['VencA_L']+self.pdk['M1']['Width']//2, check=False)))

        self.fin = self.addGen( Wire( 'fin', 'Fin', 'h',
                                      clg=UncoloredCenterLineGrid( pitch= self.pdk['Fin']['Pitch'], width= self.pdk['Fin']['Width'], offset= self.pdk['Fin']['Offset']),
                                      spg=SingleGrid( offset=0, pitch=unitCellLength)))

        stoppoint = ((self.gateDummy-1)* self.pdk['Poly']['Pitch'] +  self.pdk['Poly']['Offset'])*(1-self.shared_diff)
        self.active = self.addGen( Wire( 'active', 'Active', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=activePitch, width=activeWidth, offset=activeOffset),
                                         spg=EnclosureGrid( pitch=unitCellLength, offset=0, stoppoint=stoppoint, check=True)))

        self.RVT = self.addGen( Wire( 'RVT', 'Rvt', 'h',
                                      clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                      spg=EnclosureGrid( pitch=unitCellLength, offset=0, stoppoint=stoppoint, check=True)))

        self.LVT = self.addGen( Wire( 'LVT', 'Lvt', 'h',
                                      clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                      spg=EnclosureGrid( pitch=unitCellLength, offset=0, stoppoint=stoppoint, check=True)))

        self.HVT = self.addGen( Wire( 'HVT', 'Hvt', 'h',
                                      clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                      spg=EnclosureGrid( pitch=unitCellLength, offset=0, stoppoint=stoppoint, check=True)))

        self.SLVT = self.addGen( Wire( 'SLVT', 'Slvt', 'h',
                                      clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                      spg=EnclosureGrid( pitch=unitCellLength, offset=0, stoppoint=stoppoint, check=True)))

        offset = self.gateDummy*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Offset'] - self.pdk['Poly']['Pitch']//2
        stoppoint = self.gateDummy*self.pdk['Poly']['Pitch'] + self.pdk['Poly']['Offset']-self.pdk['Pc']['PcExt']-self.pdk['Poly']['Width']//2

        # Calculate Pc width to ensure proper V0 enclosure
        # Pc width must be >= V0 height + 2 * V0 enclosure by Poly requirement
        v0_width_y = self.pdk['V0']['WidthY']  # 90nm
        v0_enc_poly = self.pdk['V0']['VencP_L']  # 40nm
        pc_width_for_v0 = v0_width_y + 2 * v0_enc_poly  # 170nm
        pc_width = max(self.pdk['Pc'].get('PcWidth', 100), pc_width_for_v0)

        self.pc = self.addGen( Wire( 'pc', 'Pc', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=self.pdk['M2']['Pitch'], width=pc_width, offset=self.pdk['M2']['Pitch']),
                                         spg=EnclosureGrid( pitch=unitCellLength, offset=offset*self.shared_diff, stoppoint=stoppoint-offset*self.shared_diff, check=True)))

        self.nselect = self.addGen( Region( 'nselect', 'Nselect',
                                            v_grid=UncoloredCenterLineGrid( offset= 0, pitch= self.pdk['M3']['Pitch'], width= self.pdk['M3']['Width']),
                                            h_grid=self.fin.clg))
        self.pselect = self.addGen( Region( 'pselect', 'Pselect',
                                            v_grid=UncoloredCenterLineGrid( offset= 0, pitch= self.pdk['M3']['Pitch'], width= self.pdk['M3']['Width']),
                                            h_grid=self.fin.clg))
        self.nwell = self.addGen( Region( 'nwell', 'Nwell',
                                            v_grid=UncoloredCenterLineGrid( offset= 0, pitch= self.pdk['M3']['Pitch'], width= self.pdk['M3']['Width']),
                                            h_grid=self.fin.clg))

        # Add outline region to define cell boundary
        self.outline = self.addGen( Region( 'outline', 'Outline',
                                            v_grid=UncoloredCenterLineGrid( offset= 0, pitch= self.pdk['M3']['Pitch'], width= self.pdk['M3']['Width']),
                                            h_grid=self.fin.clg))

        self.active_diff = self.addGen( Wire( 'active_diff', 'Active', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=activePitch, width=activeWidth, offset=activeOffset),
                                         spg=SingleGrid( pitch=self.pdk['Poly']['Pitch'], offset=(self.gateDummy-1)*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Pitch']//2)))

        self.RVT_diff = self.addGen( Wire( 'RVT_diff', 'Rvt', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                         spg=SingleGrid( pitch=self.pdk['Poly']['Pitch'], offset=(self.gateDummy-1)*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Pitch']//2)))

        self.LVT_diff = self.addGen( Wire( 'LVT_diff', 'Lvt', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                         spg=SingleGrid( pitch=self.pdk['Poly']['Pitch'], offset=(self.gateDummy-1)*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Pitch']//2)))

        self.HVT_diff = self.addGen( Wire( 'HVT_diff', 'Hvt', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                         spg=SingleGrid( pitch=self.pdk['Poly']['Pitch'], offset=(self.gateDummy-1)*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Pitch']//2)))

        self.SLVT_diff = self.addGen( Wire( 'SLVT_diff', 'Slvt', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=activePitch, width=RVTWidth, offset=activeOffset),
                                         spg=SingleGrid( pitch=self.pdk['Poly']['Pitch'], offset=(self.gateDummy-1)*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Pitch']//2)))

        stoppoint = unitCellLength//2-self.pdk['Active']['activebWidth_H']//2
        offset_active_body = (self.lFin//2)*self.pdk['Fin']['Pitch']+self.unitCellHeight-self.pdk['Fin']['Pitch']//2
        self.activeb = self.addGen( Wire( 'activeb', 'Active', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=self.pdk['M2']['Pitch'], width=self.pdk['Active']['activebWidth'], offset=0),
                                         spg=EnclosureGrid( pitch=unitCellLength, offset=0, stoppoint=stoppoint, check=True)))

        self.activeb_diff = self.addGen( Wire( 'activeb_diff', 'Active', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=self.pdk['M2']['Pitch'], width=self.pdk['Active']['activebWidth'], offset=0),
                                         spg=SingleGrid( pitch=self.pdk['Poly']['Pitch'], offset=(self.gateDummy-1)*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Pitch']//2)))

        self.pb_diff = self.addGen( Wire( 'pb_diff', 'Pb', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=self.pdk['M2']['Pitch'], width=self.pdk['Pb']['pbWidth'], offset=0),
                                         spg=SingleGrid( pitch=self.pdk['Poly']['Pitch'], offset=(self.gateDummy-1)*self.pdk['Poly']['Pitch']+self.pdk['Poly']['Pitch']//2)))

        stoppoint = unitCellLength//2-self.pdk['Pb']['pbWidth_H']//2
        self.pb = self.addGen( Wire( 'pb', 'Pb', 'h',
                                         clg=UncoloredCenterLineGrid( pitch=self.pdk['M2']['Pitch'], width=self.pdk['Pb']['pbWidth'], offset=0),
                                         spg=EnclosureGrid( pitch=unitCellLength, offset=0, stoppoint=stoppoint, check=True)))


        self.v1_x = self.addGen( Via( 'v1_x', 'V1',
                                    h_clg=self.m2_updated.clg,
                                    v_clg=self.m1_updated.clg,
                                    WidthX=self.pdk['V1']['WidthX'],
                                    WidthY=self.pdk['V1']['WidthY']))

        self.va = self.addGen( Via( 'va', 'V0',
                                    h_clg=self.m2.clg,
                                    v_clg=self.m1_updated.clg,
                                    WidthX=self.pdk['V0']['WidthX'],
                                    WidthY=self.pdk['V0']['WidthY']))

        self.v0 = self.addGen( Via( 'v0', 'V0',
                                    h_clg=CenterLineGrid(),
                                    v_clg=self.m1_updated.clg,
                                    WidthX=self.pdk['V0']['WidthX'],
                                    WidthY=self.pdk['V0']['WidthY']))

        self.v0.h_clg.addCenterLine( 0,                 self.pdk['V0']['WidthY'], False)
        v0pitch = self.pdk['V0']['WidthY'] + self.pdk['V0']['SpaceY']
        v0Offset = activeOffset - activeWidth//2 + self.pdk['V0']['VencA_L'] + self.pdk['V0']['WidthY']//2
        v0_number = floor((activeWidth -  2*self.pdk['V0']['VencA_L'] - self.pdk['V0']['WidthY']) / v0pitch) + 1
        v0_number = max(v0_number, 1)
        assert v0_number > 0, "V0 can not be placed in the active region"
        for i in range(v0_number):
            self.v0.h_clg.addCenterLine(i*v0pitch+v0Offset,    self.pdk['V0']['WidthY'], True)
        self.v0.h_clg.addCenterLine( self.unitCellHeight,    self.pdk['V0']['WidthY'], False)
        info = self.pdk['V0']

    def _addMOS( self, x, y, x_cells,  vt_type, name='M1', reflect=False, **parameters):

        fullname = f'{name}_X{x}_Y{y}'
        self.subinsts[fullname].parameters.update(parameters)

        def _connect_diffusion(i, pin):
            # Adjust M1 length to cover V0 but avoid shorting with Gate M1
            self.addWire( self.m1_updated, f'{fullname}:{pin}', i, (grid_y0, -1), (grid_y1+1, -1), netType='pin')
            for j in range(1,self.v0.h_clg.n): ## self.v0.h_clg.n??
                self.addVia( self.v0, f'{fullname}:{pin}', i, (y, j))
            self._xpins[name][pin].append(i)

        # Draw FEOL Layers
        if self.shared_diff == 0:
            self.addWire( self.active, None, y, (x,1), (x+1,-1))
        elif self.shared_diff == 1 and x == x_cells-1:
            self.addWire( self.active_diff, None, y, 0, self.gate*x_cells+1)
        else:
            pass

        for i in range(self.gate):
            self.addWire( self.pl, None, i+self.gatesPerUnitCell*x+self.gateDummy,   (y,1), (y+1,-1))

        def _addRVT(x, y, x_cells):
            if self.shared_diff == 0:
                self.addWire( self.RVT,  None, y,          (x, 1), (x+1, -1))
            elif self.shared_diff == 1 and x == x_cells-1:
                self.addWire( self.RVT_diff,  None, y, 0, self.gate*x_cells+1)
            else:
                pass

        def _addLVT(x, y, x_cells):
            if self.shared_diff == 0:
                self.addWire( self.LVT,  None, y,          (x, 1), (x+1, -1))
            elif self.shared_diff == 1 and x == x_cells-1:
                self.addWire( self.LVT_diff,  None, y, 0, self.gate*x_cells+1)
            else:
                pass

        def _addHVT(x, y, x_cells):
            if self.shared_diff == 0:
                self.addWire( self.HVT,  None, y,          (x, 1), (x+1, -1))
            elif self.shared_diff == 1 and x == x_cells-1:
                self.addWire( self.HVT_diff,  None, y, 0, self.gate*x_cells+1)
            else:
                pass
        if vt_type == 'RVT':
            _addRVT(x, y, x_cells)
        elif vt_type == 'LVT':
            _addLVT(x, y, x_cells)
        elif vt_type == 'HVT':
            _addHVT(x, y, x_cells)
        else:
            print("This VT type not supported")
            exit()

        # Source, Drain, Gate Connections
        grid_y0 = y*self.m2PerUnitCell + 1
        grid_y1 = (y+1)*self.m2PerUnitCell-5
        gate_x = self.gateDummy*self.shared_diff + x * self.gatesPerUnitCell + self.gatesPerUnitCell // 2
        # Connect Gate (gate_x)
        self.addWire( self.m1_updated, f'{fullname}:G', gate_x , (grid_y1+2, -1), (grid_y1+3, 1), netType='pin')
        self.addWire( self.pc, None, grid_y1+1, (x,1), (x+1,-1))
        self.addVia( self.va, f'{fullname}:G', gate_x, grid_y1+2)
        self._xpins[name]['G'].append(gate_x)

        # Connect Source & Drain
        (center_terminal, side_terminal) = ('S', 'D') if self.gate%4 == 0 else ('D', 'S')
        center_terminal = 'D' if self.stack > 1 else center_terminal # for stacked devices
        _connect_diffusion(gate_x, center_terminal)

        for x_terminal in range(1,1+self.gate//2):
            terminal = center_terminal if x_terminal%2 == 0 else side_terminal #for fingers
            if self.stack > 1 and x_terminal != self.gate//2: continue
            terminal = 'S' if self.stack > 1 else terminal
            if reflect:
                _connect_diffusion(gate_x - x_terminal, terminal)
                _connect_diffusion(gate_x + x_terminal, terminal)
            else:
                _connect_diffusion(gate_x - x_terminal, terminal)
                _connect_diffusion(gate_x + x_terminal, terminal)

    def _connectDevicePins(self, y, y_cells, connections):
        center_track = y * self.m2PerUnitCell + self.m2PerUnitCell // 2 # Used for m1 extension
        grid_y1 = (y+1)*self.m2PerUnitCell-5
        gate_track = 0
        diff_track = 1
        for (net, conn) in connections.items():
            contactsx = {(track, pin) for inst, pins in self._xpins.items()
                              for pin, m1tracks in pins.items()
                              for track in m1tracks if (inst, pin) in conn}
            pins = {x[1] for x in contactsx}
            print(f"[DEBUG] Net={net}, Pins={pins}, All contactsx={contactsx}")
            for pin in pins:
                contacts = {x[0] for x in contactsx if x[1]==pin}
                print(f"[DEBUG] Pin={pin}, contacts count={len(contacts)}, contacts={sorted(contacts)}")

                # For Gate pin, only connect center M1 track to reduce M2 width
                if pin == 'G' and len(contacts) > 1:
                    print(f"[DEBUG] Gate pin: original contacts count = {len(contacts)}, contacts = {sorted(contacts)}")
                    sorted_contacts = sorted(list(contacts))
                    mid_index = len(sorted_contacts) // 2
                    contacts = {sorted_contacts[mid_index]}
                    print(f"[DEBUG] Gate pin: filtered to 1 contact, contacts = {contacts}")

                for j in range(self.minvias):
                    if pin == 'G':
                        current_track = grid_y1 + 2 + gate_track
                        gate_track = gate_track+1
                    elif pin == 'B':
                        if y != y_cells-1:continue
                        current_track = y_cells * self.m2PerUnitCell + self.lFin//4
                    else:
                        current_track = y * self.m2PerUnitCell + len(connections) * j + diff_track
                        diff_track = diff_track + 1
                    print(f"[DEBUG] addWireAndViaSet: net={net}, pin={pin}, current_track={current_track}, contacts={contacts}")
                    self.addWireAndViaSet(net, self.m2_updated, self.v1_x, current_track, contacts)
                    self._nets[net][current_track] = contacts
                # Extend m1 if needed. TODO: Should we draw longer M1s to begin with?
                #direction = 1 if current_track > center_track else -1
                #for i in contacts:
                #    self.addWire( self.m1, net, None, i, (center_track, -1 * direction), (current_track, direction))


    def _connectNets(self, x_cells, y_cells):
        def _get_wire_terminators(intersecting_tracks):
            minx, maxx = min(intersecting_tracks), max(intersecting_tracks)
            # BEGIN: Quick & dirty MinL DRC error fix.
            minL = 2 # TODO: Make this MinL dependent
            L = maxx - minx
            if center_track - minL // 2 <= minx and maxx <= center_track + minL // 2:
                minx, maxx = (center_track - minL // 2, center_track + minL // 2)
            elif L < minL:
                minx, maxx = (minx - (minL - L), maxx) if minx >= center_track else (minx, maxx + (minL - L))
            # END: Quick & dirty MinL DRC error fix.
            return (minx, maxx)
        M1_tracks = x_cells*self.gatesPerUnitCell+2*self.gateDummy*self.shared_diff
        M3_tracks = ceil(M1_tracks*self.pdk['M1']['Pitch']/self.pdk['M3']['Pitch'])
        center_track = (M3_tracks - len(self._nets) * self.minvias)// 2
        m3start = center_track
        #m3start = self.shared_diff+(x_cells * self.gatesPerUnitCell - len(self._nets) * self.minvias) // 2
        for track, (net, conn) in enumerate(self._nets.items(), start=1):
            for j in range(self.minvias):
                current_track = m3start + len(self._nets) * j + track
                contacts = conn.keys()
                if len(contacts) == 1: # Create m2 terminal
                    i = next(iter(contacts))
                    # For Gate pin, don't expand M2 width with MinL rule
                    if net == 'G':
                        m1_tracks = list(conn[i])
                        minx, maxx = min(m1_tracks), max(m1_tracks)
                        print(f"[DEBUG] Gate M2 pin: M2_track={i}, M1_tracks={m1_tracks}, minx={minx}, maxx={maxx}")
                    else:
                        minx, maxx = _get_wire_terminators(conn[i])
                    self.addWire(self.m2_updated, net, i, (minx, -1), (maxx, 1), netType = 'pin')
                else: # create m3 terminal(s)
                    self.addWireAndViaSet(net, self.m3, self.v2, current_track, contacts, netType = 'pin')
                    if max(contacts) - min(contacts) < 2:
                        minL_m3 = 2 ## for M3 min length
                        maxy = min(contacts) + minL_m3
                        miny = min(contacts)
                        self.addWire(self.m3, net, current_track, (miny, -1), (maxy, 1), netType = 'pin')
                    else:
                        pass
                    # Extend m2 if needed. TODO: What to do if we go beyond cell boundary?
                    for i, locs in conn.items():
                        minx, maxx = _get_wire_terminators([*locs, current_track])
                        if self.pdk['M3']['Pitch'] >= self.pdk['M1']['Pitch']:
                            self.addWire(self.m2, net, i, (minx, -1), (maxx, 1))
                        else:
                            self.addWire( self.m2_updated, net, i, (0, 1), (M1_tracks, -1))

    def _addBodyContact(self, x, y, x_cells, yloc=None, name='M1'):
        fullname = f'{name}_X{x}_Y{y}'
        if yloc is not None:
            y = yloc
        h = self.m2PerUnitCell
        gu = self.gatesPerUnitCell
        body_v0_track = (self.lFin*self.pdk['Fin']['Pitch'])//(2*self.pdk['M2']['Pitch'])
        gate_x = self.gateDummy*self.shared_diff + x*gu + gu // 2
        self._xpins[name]['B'].append(gate_x)
        if self.shared_diff == 0:
            self.addWire( self.activeb, None, (y+1)*h + body_v0_track, (x,1), (x+1,-1))
            self.addWire( self.pb, None, (y+1)*h + body_v0_track, (x,1), (x+1,-1))
        else:
            self.addWire( self.activeb_diff, None, (y+1)*h + body_v0_track, 0, self.gate*x_cells+1)
            self.addWire( self.pb_diff, None, (y+1)*h + body_v0_track, (x,1), (x+1,-1))
        self.addWire( self.m1_updated, f'{fullname}:B', gate_x, ((y+1)*h + body_v0_track-1, -1), ((y+1)*h + body_v0_track+1, 1), netType='pin')
        self.addVia( self.va, f'{fullname}:B', gate_x, (y+1)*h + body_v0_track)

    def _addMOSArray( self, x_cells, y_cells, pattern, vt_type, connections, minvias = 1, **parameters):
        if minvias * len(connections) > self.m2PerUnitCell - 1:
            self.minvias = (self.m2PerUnitCell - 1) // len(connections)
            logger.warning( f"Using minvias = {self.minvias}. Cannot route {len(connections)} signals using minvias = {minvias} (max m2 / unit cell = {self.m2PerUnitCell})" )
        else:
            self.minvias = minvias
        names = ['M1'] if pattern == 0 else ['M1', 'M2']
        names = sorted({c[0] for mc in connections.values() for c in mc})
        self._nets = collections.defaultdict(lambda: collections.defaultdict(list)) # net:m2track:m1contacts (Updated by self._connectDevicePins)
        ### Needs to be generalized
        if len(parameters) > 2:
            device_name_all = [*parameters.keys()]
            if int(parameters[device_name_all[0]]["NFIN"])*int(parameters[device_name_all[0]]["NF"])*int(parameters[device_name_all[0]]["M"]) != int(parameters[device_name_all[1]]["NFIN"])*int(parameters[device_name_all[1]]["NF"])*int(parameters[device_name_all[1]]["M"]):
                pattern=3
                if int(parameters[device_name_all[0]]["NF"])*int(parameters[device_name_all[0]]["M"]) > int(parameters[device_name_all[1]]["NF"])*int(parameters[device_name_all[1]]["M"]):
                    x_left = x_cells//2 - (int(parameters[device_name_all[1]]["NF"])*int(parameters[device_name_all[1]]["M"]))//2
                    x_right = x_cells//2 + (int(parameters[device_name_all[1]]["NF"])*int(parameters[device_name_all[1]]["M"]))//2
                else:
                    x_left = x_cells//2 - (int(parameters[device_name_all[0]]["NF"])*int(parameters[device_name_all[0]]["M"]))//2
                    x_right = x_cells//2 + (int(parameters[device_name_all[0]]["NF"])*int(parameters[device_name_all[0]]["M"]))//2
         ##########################

        for y in range(y_cells):
            self._xpins = collections.defaultdict(lambda: collections.defaultdict(list)) # inst:pin:m1tracks (Updated by self._addMOS)

            for x in range(x_cells):

                if self.exact_patterns: # Exact pattern from user
                    row_pattern = self.exact_patterns[y][x]
                    names_mapping = list(string.ascii_uppercase)
                    names_updated = {}
                    for i in range(len(names)):
                        names_updated[names_mapping[i]] = names[i]
                        names_updated[names_mapping[i].lower()] = names[i]
                    reflect = row_pattern.islower()
                    self._addMOS(x, y, x_cells, vt_type, names_updated[row_pattern],  False, **parameters)
                    # Removed body contact for basic MOS primitive
                    # if self.bodyswitch==1:self._addBodyContact(x, y, x_cells, y_cells - 1, names_updated[row_pattern])
                elif pattern == 0: # None (single transistor)
                    # TODO: Not sure this works without dummies. Currently:
                    # A A A A A A
                    self._addMOS(x, y, x_cells, vt_type, names[0], False, **parameters)
                    # Removed body contact for basic MOS primitive
                    # if self.bodyswitch==1:self._addBodyContact(x, y, x_cells, y_cells - 1, names[0])
                elif pattern == 1: # CC
                    # TODO: Think this can be improved. Currently:
                    # A B B A A' B' B' A'
                    # B A A B B' A' A' B'
                    # A B B A A' B' B' A'
                    self._addMOS(x, y, x_cells, vt_type, names[((x // 2) % 2 + x % 2 + (y % 2)) % 2], x >= x_cells // 2,  **parameters)
                    # Removed body contact for basic MOS primitive
                    # if self.bodyswitch==1:self._addBodyContact(x, y, x_cells, y_cells - 1, names[((x // 2) % 2 + x % 2 + (y % 2)) % 2])
                elif pattern == 2: # interdigitated
                    # TODO: Evaluate if this is truly interdigitated. Currently:
                    # A B A B A B
                    # B A B A B A
                    # A B A B A B
                    self._addMOS(x, y, x_cells, vt_type, names[((x % 2) + (y % 2)) % 2], False,  **parameters)
                    # Removed body contact for basic MOS primitive
                    # if self.bodyswitch==1:self._addBodyContact(x, y, x_cells, y_cells - 1, names[((x % 2) + (y % 2)) % 2])
                elif pattern == 3: # CurrentMirror
                    # TODO: Evaluate if this needs to change. Currently:
                    # B B B A A B B B
                    # B B B A A B B B
                    self._addMOS(x, y, x_cells, vt_type, names[0 if x_left <= x < x_right else 1], False,  **parameters)
                    # Removed body contact for basic MOS primitive
                    # if self.bodyswitch==1:self._addBodyContact(x, y, x_cells, y_cells - 1, names[0 if x_left <= x < x_right else 1])
                elif pattern == 4:  # non common centroid
                    # TODO: Evaluate if this is truly interdigitated. Currently:
                    # A A A B B B
                    # A A A B B B
                    # A A A B B B
                    self._addMOS(x, y, x_cells, vt_type, names[0 if x < (x_cells//2) else 1], False, **parameters)
                    # Removed body contact for basic MOS primitive
                    # if self.bodyswitch == 1:
                    #     self._addBodyContact(x, y, x_cells, y_cells - 1, names[0 if x < (x_cells//2) else 1])
                else:
                    assert False, "Unknown pattern"
            self._connectDevicePins(y, y_cells, connections)
        self._connectNets(x_cells, y_cells)

    def addNMOSArray( self, x_cells, y_cells, pattern, vt_type, connections, **parameters):
        from math import ceil, floor

        self._addMOSArray(x_cells, y_cells, pattern, vt_type, connections, **parameters)

        #####   Cell Boundary and Nselect Placement - Optimized X and Y   #####
        M3_tracks_end = ceil((x_cells*self.gatesPerUnitCell+2*self.gateDummy*self.shared_diff)*self.pdk['M1']['Pitch']/self.pdk['M3']['Pitch'])
        M3_tracks_start = ceil(self.pdk['M1']['Pitch']/self.pdk['M3']['Pitch'])

        # Add outline region to fix cell boundary at full M3 tracks range
        self.addRegion( self.outline, None, -M3_tracks_start, 0, M3_tracks_end, y_cells * self.finsPerUnitCell)

        # Y direction: optimize based on Active position + DRC minimum enclosure (120nm)
        fin_pitch = self.pdk['Fin']['Pitch']
        m3_pitch = self.pdk['M3']['Pitch']

        # Calculate Active boundaries in physical coordinates
        active_y_bottom = self.activeOffset - self.activeWidth//2
        active_y_top = self.activeOffset + self.activeWidth//2

        # Add DRC enclosure (120nm) and convert to fin grid for Y direction
        # Use round() to get closest grid point for minimum enclosure
        np_y_start_phys = active_y_bottom - self.NselectEnclosure  # 730nm
        np_y_end_phys = active_y_top + self.NselectEnclosure  # 1970nm

        np_y_bottom = int(round(np_y_start_phys / fin_pitch))  # 730/100=7.3 → 7
        np_y_top = int(round(np_y_end_phys / fin_pitch))  # 1970/100=19.7 → 20

        # Verify minimum enclosure is met (safety check)
        if (np_y_bottom * fin_pitch) > np_y_start_phys:
            np_y_bottom -= 1  # Expand downward if needed
        if (np_y_top * fin_pitch) < np_y_end_phys:
            np_y_top += 1  # Expand upward if needed

        # Ensure NP stays within cell boundaries
        np_y_bottom = max(0, np_y_bottom)
        np_y_top = min(y_cells * self.finsPerUnitCell, np_y_top)

        # X direction: optimize NP enclosure of Active with DRC minimum (120nm)
        # Active X position is complex to predict due to EnclosureGrid wire placement
        # Use empirically observed values and calculate symmetric enclosure

        # For NF=2 NMOS: Active empirically observed at ~1100-1940nm (840nm width)
        # Calculate Active width based on gate structure
        poly_pitch = self.pdk['Poly']['Pitch']
        active_gates = self.gate if self.shared_diff else (self.gatesPerUnitCell - 2 * self.gateDummy)

        # Active width: approximately 1.5 * poly_pitch for 2-gate device
        # Empirical observation: 840nm for NF=2 (2 gates × 280nm × 1.5 = 840nm)
        active_x_width = int(active_gates * poly_pitch * 1.5)

        # Active center is offset from cell center due to dummy gate structure
        # Cell width from M3 grid
        cell_width_physical = (M3_tracks_end + M3_tracks_start) * m3_pitch

        # Empirically, Active center is at approximately 54.3% of cell width
        # For 2800nm cell: Active center ≈ 1520nm = (1100+1940)/2
        active_center_ratio = 0.5429  # 1520/2800
        active_x_center = int(cell_width_physical * active_center_ratio)

        # Calculate Active boundaries from center and width
        active_x_start = active_x_center - active_x_width // 2
        active_x_end = active_x_center + active_x_width // 2

        # NP boundaries with DRC enclosure (120nm)
        np_x_start_phys = active_x_start - self.NselectEnclosure
        np_x_end_phys = active_x_end + self.NselectEnclosure

        # Convert physical coordinates to M3 grid indices
        # Use floor for start (expand leftward) and ceil for end (expand rightward)
        # This ensures enclosure is always >= 120nm
        np_m3_start_abs = int(floor(np_x_start_phys / m3_pitch))
        np_m3_end_abs = int(ceil(np_x_end_phys / m3_pitch))

        # Convert absolute grid to relative grid indices
        np_m3_start_index = np_m3_start_abs - M3_tracks_start
        np_m3_end_index = np_m3_end_abs - M3_tracks_start

        # Clamp to valid cell range
        np_m3_start = max(-M3_tracks_start, min(np_m3_start_index, M3_tracks_end))
        np_m3_end = max(-M3_tracks_start, min(np_m3_end_index, M3_tracks_end))

        self.addRegion( self.nselect, None, np_m3_start, np_y_bottom, np_m3_end, np_y_top)
        # Removed Pselect for NMOS - not required for basic MOS primitive
        # if self.bodyswitch==1:self.addRegion( self.pselect, None, -M3_tracks_start, y_cells* self.finsPerUnitCell, M3_tracks_end, y_cells* self.finsPerUnitCell+self.bodyswitch*self.lFin)

    def addPMOSArray( self, x_cells, y_cells, pattern, vt_type, connections, **parameters):

        self._addMOSArray(x_cells, y_cells, pattern, vt_type, connections, **parameters)

        #####   Pselect and Nwell Placement   #####
        M3_tracks_end = ceil((x_cells*self.gatesPerUnitCell+2*self.gateDummy*self.shared_diff)*self.pdk['M1']['Pitch']/self.pdk['M3']['Pitch'])
        M3_tracks_start = ceil(self.pdk['M1']['Pitch']/self.pdk['M3']['Pitch'])

        self.addRegion( self.pselect, None, -M3_tracks_start, 0, M3_tracks_end, y_cells* self.finsPerUnitCell)
        if self.bodyswitch==1:self.addRegion( self.nselect, None, -M3_tracks_start, y_cells* self.finsPerUnitCell, M3_tracks_end, y_cells* self.finsPerUnitCell+self.bodyswitch*self.lFin)
        self.addRegion( self.nwell, None, -M3_tracks_start, 0, M3_tracks_end, y_cells* self.finsPerUnitCell+self.bodyswitch*self.lFin)


