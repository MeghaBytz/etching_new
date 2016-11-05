function vdot=new_oxys_diss_iz(t,v)
global ee me MO MO2 ng Tg Ti l_p gammaO gammaO2m Efactor_O2
global Efactor_O EnergyO2 sigO2 EnergyO sigO
global Krec Krec2 Krec3 Krec4 Kdet Kch Rlambda hl0 Rrec alphabar
global pabs R area volume QtorrLit Qmolec Kpump scat_Xsec
vdot=zeros(7,1);
myV = zeros(7,1);
nO2=v(1);
nO2plusbar=v(2);
nOplusbar=v(3);
nOminusbar=v(4);
nObar=v(5);
nO2mbar=v(6);
pe=v(7);
ng=nO2+nObar+nO2mbar;
ne0=nO2plusbar+nOplusbar-nOminusbar;
Te=pe/(1.5*ne0);
Tplus=Ti;
Tminus=Ti;
gamma_plus=Te/Tplus;
gamma_minus=Te/Tminus;
% ratio for density weighting 
nOplus_ratio=nOplusbar/(nO2plusbar+nOplusbar);
% calculating lambda
% old way --> lambda=1/(330*p*1e-3)*1e-2; % in m
lambda=1/(ng*scat_Xsec); % lambda in m
% caculating values for O neutral & O2m
vbarO=sqrt(8*ee*Ti/(pi*MO)); % average thermal velocity of O neutral
vbarO2m=sqrt(8*ee*Ti/(pi*MO2)); % average thermal velocity of O2m
DO=ee*Tg*lambda/vbarO/MO; % Diffusion coefficient of O neutral
DO2m=ee*Tg*lambda/vbarO2m/MO2; % Diffusion coefficient of O2m
dO=sqrt(4*DO*l_p*(2-gammaO)/vbarO/gammaO + l_p^2);
dO2m=sqrt(4*DO2m*l_p*(2-gammaO2m)/vbarO2m/gammaO2m + l_p^2);
hAO=1/(1 + l_p*vbarO*gammaO/4/DO/(2-gammaO));
hAO2m=1/(1 + l_p*vbarO2m*gammaO2m/4/DO2m/(2-gammaO2m));
% wall recombination loss rate for O neutral
KO=hAO*vbarO*2*gammaO/(2-gammaO)*area/volume/4;
% wall recombination loss rate for O2m
KO2m=hAO2m*vbarO2m*2*gammaO2m/(2-gammaO2m)*area/volume/4;
vol_O=volume*(1 - l_p^2/(3*dO^2))*(1 - (2/3)*l_p^3/(R*dO^2)...
 + l_p^4/(6*R^2*dO^2));
vol_O2m=volume*(1 - l_p^2/(3*dO2m^2))*(1 - (2/3)*l_p^3/(R*dO2m^2)...
 + l_p^4/(6*R^2*dO2m^2));
%finding alpha0
eta=2*Tplus/(Tplus+Tminus);
hl0=0.86/sqrt(3+eta*l_p/lambda);
Rlambda=sqrt(2*pi/gamma_plus)*lambda/l_p/eta;
uB_O2=sqrt(ee*Te/MO2);
uB_O=sqrt(ee*Te/MO);
% density-weighted Bohm velocity
uB_dw=uB_O2*(1-nOplus_ratio)+uB_O*nOplus_ratio;
% density-weighted Recombination rate coefficient
Krec_dw=(Krec+Krec2)*(1-nOplus_ratio)+Krec3*nOplus_ratio;
Rrec=Krec_dw*ne0*l_p/uB_dw;
alphabar=nOminusbar/ne0;
alpha0=fzero('alpha0find',[1e-2 1e2]);
volminus=vol_minus(alpha0);
%calculating peak(center) densities
nOminus=nOminusbar*volume/volminus;
nplus=nOminus+ne0;
nOplus=nplus*nOplus_ratio;
nO2plus=nplus-nOplus;
nO=nObar*volume/vol_O;
nO2m=nO2mbar*volume/vol_O2m;
% calculation of hl factors with 3 models
% density weighted positive ion mass
Mplus_dw=MO2*(1-nOplus_ratio)+MO*nOplus_ratio;
nstar=15/56*sqrt(8*ee*Tplus/pi/Mplus_dw)*(eta^2)/(Krec_dw*lambda);
hpar2=hl0*1/(1+alpha0);
hpar1=1/(gamma_minus^0.5 + (gamma_plus^0.5)*(eta*l_p/sqrt(2*pi)/lambda))*(alpha0/(1+alpha0));
hflat1=1/(gamma_minus^0.5 + (gamma_plus^0.5)*(nstar^0.5)/(nOminus^0.5));
hl=hpar2 + hpar1 + hflat1; %overall hl factor by "linear ansatz"
% size of EN core (lminus & rminus)
l_minus=lminus(alpha0);
r_minus=R - l_p + l_minus;
% Effective volume for recombination losses
vol_rec=(2*pi*l_minus/(1+alpha0))*(8/15*alpha0*(r_minus^2 - 14/15*r_minus*l_minus + 4/15*l_minus^2) + 2/3*(r_minus^2 - 2/3*r_minus*l_minus + 1/6*l_minus^2));
% Surface loss of electron-ion pair
Kion=hl*uB_dw*area/volume;
% other reaction coefficients
%Kiz1=2.34E-15*Te^(1.03)*exp(-12.29/Te);
Kei=2.2E-14*Te^(-0.5);
Katt=1.07E-15*Te^(-1.391)*exp(-6.26/Te);
Kiz2=9E-15*Te^(0.7)*exp(-13.6/Te);
Kdiss=6.09*6.86E-15*exp(-6.29/Te);
Kiz3=7.1E-17*Te^0.5*exp(-17/Te);
Kiz4=1.88E-16*Te^(1.699)*exp(-16.81/Te);
Kex=1.37E-15*exp(-2.14/Te);
Kizm=2.34E-15*Te^(1.03)*exp(-11.31/Te);
Kattm=4.19E-15*Te^(-1.376)*exp(-5.19/Te);
Kdeex=2.06E-15*exp(-1.163/Te);
Kdism=6.09*6.86E-15*exp(-5.31/Te);
% calculating Ec
Ec_O=o_ec(Te);
Ec_O2=o2_ec(Te);
Eei_O=Efactor_O*Te;
Eei_O2=Efactor_O2*Te;
% calculating reduction factor for volume of ionization
Kel=rateconstant( EnergyO2, sigO2, Te );
vbare=sqrt(8*ee*Te/(pi*me)); % average thermal velocity of electron
lambda_E=vbare/ng/sqrt(3*Kel*Kex);
vr_iz=1/(1 + 2*l_p/lambda_E);
Kiz1=vr_iz*2.34E-15*Te^(1.03)*exp(-12.29/Te); % reduced Kiz1
% differential EQ's
% for nO2
vdot(1)=Qmolec/volume + Krec*nO2plus*nOminus*vol_rec/volume...
 + Kdet*nOminusbar*nObar + Kdeex*nO2mbar*ne0...
 + Krec4*nO2mbar*nOminusbar + Kion*nO2plus + KO2m*nO2m...
 + 0.5*KO*nO - (Kiz1 + Katt + Kdiss + Kiz3 + Kiz4 + Kex)*nO2*ne0...
 - Kch*nOplusbar*nO2 - Kpump*nO2;
% for nO2plus
vdot(2)=Kiz1*nO2*ne0 + Kizm*nO2mbar*ne0 + Kch*nOplusbar*nO2...
 -(Krec + Krec2)*nO2plus*nOminus*vol_rec/volume...
 - Kei*ne0*nO2plusbar - Kion*nO2plus;
% for nOplus
vdot(3)=Kiz2*nObar*ne0 + (Kiz3 + Kiz4)*nO2*ne0...
 - Krec3*nOplus*nOminus*vol_rec/volume - Kch*nOplusbar*nO2...
 - Kion*nOplus;
% for nOminus
vdot(4)=(Katt+Kiz3)*nO2*ne0 + Kattm*nO2mbar*ne0...
 - ((Krec+Krec2)*nO2plus*nOminus+Krec3*nOplus*nOminus)...
 *vol_rec/volume - Kdet*nOminusbar*nObar - Krec4*nOminusbar*nO2mbar;
% for nO
vdot(5)=2*Kei*ne0*nO2plusbar + (2*Kdiss+Katt+Kiz4)*ne0*nO2...
 + (Krec+3*Krec2)*nO2plus*nOminus*vol_rec/volume...
 + 2*Krec3*nOplus*nOminus*vol_rec/volume + Kch*nOplusbar*nO2...
 + (Kattm+2*Kdism)*nO2mbar*ne0 + Krec4*nOminusbar*nO2mbar...
 + Kion*nOplus - Kiz2*nObar*ne0 - Kdet*nOminusbar*nObar...
 - KO*nO - Kpump*nObar;
% for nO2m
vdot(6)=Kex*nO2*ne0 - (Kizm+Kattm+Kdeex+Kdism)*nO2mbar*ne0...
 - Krec4*nO2mbar*nOminusbar - KO2m*nO2m - Kpump*nO2mbar; 
% for power balance
vdot(7)=pabs - Ec_O2*Kiz1*nO2*ne0 - Ec_O*Kiz2*nObar*ne0...
 - Eei_O*Kion*nOplus - Eei_O2*Kion*nO2plus; 

%define my variables
k1=	Kiz1;
k2=	Katt;
k3=	Kiz4;
k4=	Kiz3;
k5=	Kdiss;
k6=	0;
k7=	Kiz2;
k8=	Kei;
k9=	Kdet;
k10=Kch;
k11=Kex;
k12=Kdeex;
k13=Kizm;
k14=Kattm;
k15=Kdism;
k16=Krec4;
k17=Krec;
k18=Krec2;
k19=Krec3;
k20=KO;
k21=KO2m;
k22=Kion;
k23=Kion;
n1 = nO2;
n2 = nObar;
n3 = nO2plusbar;
n4 = nOplusbar;
n5 = nOminusbar;
n6 = nO2mbar;
%n7 = nOmbar;
n8 = ne0;

np2 = nO;
np3 = nO2plus;
np4 = nOplus;
np5 = nOminus;
np6 = nO2m;
np9 = 1;
Vrec = vol_rec;
V = volume;
Q1 = Qmolec;
Ec_1 = Ec_O2;
Ec_2 = Ec_O;
Ew_1 = Eei_O2;
Ew_2 = Eei_O;
a = 1;
b = 2;

myV(1)=-k1 * n1 * n8 + -k2 * n1 * n8 + -k3 * n1 * n8 + -k4 * n1 * n8 + -k5 * n1 * n8 + -k6 * n1 * n8 + k9 * n2 * n5 + -k10 * n1 * n4 + -k11 * n1 * n8 + k12 * n6 * n8 + k16 * n5 * n6 - Kpump*n1 + k17 * np3 * np5 * Vrec/V + k20/2 * np2 * np9 + k21 * np6 * np9 + k22 * np3 * np9 + Q1/V;
myV(2)=k1 * n1 * n8 + -k8 * n3 * n8 + k10 * n1 * n4 + k13 * n6 * n8 + -k17 * np3 * np5 * Vrec/V + -k18 * np3 * np5 * Vrec/V + -k22 * np3 * np9;
myV(3)=k3 * n1 * n8 + k4 * n1 * n8 + k7 * n2 * n8 + -k10 * n1 * n4 + -k19 * np4 * np5 * Vrec/V + -k23 * np4 * np9;
myV(4)=k2 * n1 * n8 + k4 * n1 * n8 + -k9 * n2 * n5 + k14 * n6 * n8 + -k16 * n5 * n6 + -k17 * np3 * np5 * Vrec/V + -k18 * np3 * np5 * Vrec/V + -k19 * np4 * np5 * Vrec/V;
myV(5)=k2 * n1 * n8 + k3 * n1 * n8 + 2*k5 * n1 * n8 + k6 * n1 * n8 + -k7 * n2 * n8 + 2*k8 * n3 * n8 + -k9 * n2 * n5 + k10 * n1 * n4 + k14 * n6 * n8 + 2*k15 * n6 * n8 + k16 * n5 * n6 - Kpump*n2 + k17 * np3 * np5 * Vrec/V + 3*k18 * np3 * np5 * Vrec/V + 2*k19 * np4 * np5 * Vrec/V + -k20 * np2 * np9 + k23 * np4 * np9;
myV(6)=k11 * n1 * n8 + -k12 * n6 * n8 + -k13 * n6 * n8 + -k14 * n6 * n8 + -k15 * n6 * n8 + -k16 * n5 * n6 - Kpump*n6 + -k21 * np6 * np9;
%myV(7)=k6 * n1 * n8 - Kpump*n7;
myV(7)=pabs - Ec_1*k1 * n1 * n8 - Ec_2*k7 * n2 * n8 - Ew_1*k22 * np3 * np9 - Ew_2*k23 * np4 * np9;