function oEc=o_ec(Te); 
global me MO EnergyO sigO
Energy=EnergyO;
sig=sigO;
Kel=rateconstant( Energy, sig, Te );
kiz = 9e-15 * (Te^0.7) * exp(-13.6/Te);
k1D = 4.54e-15 * exp(-2.36/Te);
k1S = 7.86e-16 * exp(-4.489/Te);
k3P0 = 2.53e-15 * exp(-17.34/Te);
k5S0 = 9.67e-16 * exp(-9.97/Te);
k3S0 = 3.89e-15 * exp(-9.75/Te);
kh = 4.31e-14 * exp(-18.59/Te);
oEc = 13.61 + k1D/kiz*1.96 + k1S/kiz*4.18 + k5S0/kiz*9.14...
 + k3S0/kiz*9.51 + k3P0/kiz*15.65 + kh/kiz*12 + Kel/kiz*3*me/MO*Te; 