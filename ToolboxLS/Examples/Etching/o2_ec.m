function o2Ec=o2_ec(Te);
global me MO2 EnergyO2 sigO2
Tel=Te;
Energy=EnergyO2;
sig=sigO2;
Kel=rateconstant( Energy, sig, Tel );
kiz = 2.34e-15*(Tel^1.03)*exp(-12.29/Tel);
krot = 1.8736e-17 * exp(-2.9055/Tel);
kv1 = 2.8e-15 * exp(-3.72/Tel);
kv2 = 1.28e-15 * exp(-3.67/Tel);
ka1D = 1.37e-15 * exp(-2.14/Tel);
kb1S = 3.24e-16 * exp(-2.218/Tel);
kex1 = 1.13e-15 * exp(-3.94/Tel);
k1dis = 6.86e-15 * exp(-6.29/Tel);
k2dis = 3.4879e-14 * exp(-5.92/Tel);
k3dis = 1.443e-16 * exp(-17.25/Tel);
kex2 = 1.13e-15 * exp(-18.35/Tel);
o2Ec = 12.14 + krot/kiz*0.02 + kv1/kiz*0.19 + kv2/kiz*0.38...
 + ka1D/kiz*0.977 + kb1S/kiz*1.627 + kex1/kiz*4.5...
 + k1dis/kiz*6 + k2dis/kiz*8.4 + k3dis/kiz*9.97...
 + kex2/kiz*14.7 + Kel/kiz*3*me/MO2*Tel; 