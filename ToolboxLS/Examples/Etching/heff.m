function heffective=heff(alpha)
global Rlambda hl0
h_l0=hl0;
Rlam=Rlambda;
heffective=h_l0 + (1 - h_l0)*alpha*Rlam; 