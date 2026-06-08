.include ./models.sp

* Inverter Test Subcircuit
.subckt INVERTER_TEST vdd in out gnd
* CMOS inverter
* PMOS pull-up transistor
MP out in vdd vdd pmos_rvt W=2u L=100n
* NMOS pull-down transistor
MN out in gnd gnd nmos_rvt W=1u L=100n
.ends INVERTER_TEST
