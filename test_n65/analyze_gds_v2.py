import gdspy
import os

gds_file = 'test_n65/NMOS_TEST_0.gds'
if not os.path.exists(gds_file):
    print(f"Error: {gds_file} not found")
    exit(1)

lib = gdspy.GdsLibrary(infile=gds_file)
print(f"Library unit: {lib.unit}, precision: {lib.precision}")
print(f"Cells in library: {lib.cells.keys()}")

target_cell = None
for name in lib.cells.keys():
    if name.startswith('NMOS_S_'):
        target_cell = lib.cells[name]
        break

if target_cell:
    cell = target_cell
else:
    print("Could not find NMOS_S_ cell, using top level")
    cell = lib.top_level()[0]

print(f"Analyzing cell: {cell.name}")

poly_layer = 17
active_layer = 6
m1_layer = 31
contact_layer = 30

polys = []
actives = []
m1s = []
contacts = []

# gdspy polygons are in cell.polygons
for p in cell.polygons:
    if p.layers[0] == poly_layer:
        polys.append(p)
    elif p.layers[0] == active_layer:
        actives.append(p)
    elif p.layers[0] == m1_layer:
        m1s.append(p)
    elif p.layers[0] == contact_layer:
        contacts.append(p)

print(f"Found {len(polys)} Poly shapes")
print(f"Found {len(actives)} Active shapes")
print(f"Found {len(m1s)} M1 shapes")
print(f"Found {len(contacts)} Contact shapes")

if not actives:
    print("No Active shapes found!")
    exit()

for i, active in enumerate(actives):
    print(f"Active {i} BBox: {active.get_bounding_box()}")

# Use the first active for intersection check
active = actives[0]
active_bbox = active.get_bounding_box()
active_y_min, active_y_max = active_bbox[0][1], active_bbox[1][1]
active_x_min, active_x_max = active_bbox[0][0], active_bbox[1][0]

for i, poly in enumerate(polys):
    poly_bbox = poly.get_bounding_box()
    # Check intersection - simple bbox check first
    if (poly_bbox[0][0] < active_x_max and poly_bbox[1][0] > active_x_min and
        poly_bbox[0][1] < active_y_max and poly_bbox[1][1] > active_y_min):
        
        print(f"Poly {i} (Gate) BBox: {poly_bbox}")
        ext_top = poly_bbox[1][1] - active_y_max
        ext_bottom = active_y_min - poly_bbox[0][1]
        print(f"  Extension Top: {ext_top:.4f}")
        print(f"  Extension Bottom: {ext_bottom:.4f}")

for i, m1 in enumerate(m1s):
    m1_bbox = m1.get_bounding_box()
    print(f"M1 {i} BBox: {m1_bbox}")
    if (m1_bbox[0][0] < active_x_max and m1_bbox[1][0] > active_x_min and
        m1_bbox[0][1] < active_y_max and m1_bbox[1][1] > active_y_min):
         print(f"  M1 {i} overlaps Active")
         ext_top = m1_bbox[1][1] - active_y_max
         ext_bottom = active_y_min - m1_bbox[0][1]
         ext_right = m1_bbox[1][0] - active_x_max
         ext_left = active_x_min - m1_bbox[0][0]
         print(f"  Extension Top: {ext_top:.4f}")
         print(f"  Extension Bottom: {ext_bottom:.4f}")
         print(f"  Extension Right: {ext_right:.4f}")
         print(f"  Extension Left: {ext_left:.4f}")