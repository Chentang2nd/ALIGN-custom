.include ./models.sp

* PMOS Test Subcircuit
.subckt PMOS_TEST vdd in out gnd
* PMOS transistor with load
* M<name> <drain> <gate> <source> <bulk> <model> W=<width> L=<length>
M1 out in vdd vdd pmos_rvt W=2u L=100n
RLoad gnd out 10k
.ends PMOS_TEST
