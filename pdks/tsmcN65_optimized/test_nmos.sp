.include ./models.sp

* NMOS Test Subcircuit
.subckt NMOS_TEST vdd in out gnd
* NMOS transistor with load
* M<name> <drain> <gate> <source> <bulk> <model> W=<width> L=<length>
M1 out in gnd gnd nmos_rvt W=200n L=100n
RLoad vdd out 10k
.ends NMOS_TEST
