clear
global expConditions
%global data
global noUnknowns
global numberOfIons
global numberOfRadicals
global noExperimentalMetrics
global levelSetUnknowns
global numberOfMaterials

global ee me MO MO2 Tg Ti l_p gammaO gammaO2m
global Efactor_O2 Efactor_O EnergyO2 sigO2 EnergyO sigO
global Krec Krec2 Krec3 Krec4 Kdet Kch 
global pabs R area volume QtorrLit Qmolec Kpump scat_Xsec

noExperimentalMetrics = 2;
numberOfIons = 2;
numberOfRadicals = 2;
levelSetUnknowns = 2;
numberOfMaterials = 1;
noUnknowns = numberOfIons*2+numberOfRadicals*2+levelSetUnknowns+1;

load allSynExpConditions

expConditions = allSynExpConditions;

%make fake real parameters for novel resist material "Algernon"
A1 = [.9 2];
B1 = [3.1 5.7];
rxnProb1 = [0.4 .7];
SCl1 = [.4 .7];
lambda1 = 1;
lambda2 = .6;
epsS1 = 50;
epsD1 = 8;
noise = .02;

A2 = [.2 .1];
B2 = [1 1];
rxnProb2 = [0.4 .8];
SCl2 = [.9 .9];
epsS2 = 11;
epsD2 = 2;

material1 = [A1 B1 rxnProb1 SCl1 lambda1 lambda2 epsS1 epsD1 noise];
material2 = [A2 B2 rxnProb2 SCl2 lambda1 lambda2 epsS2 epsD2 noise];

algernonUnknowns = [material1; material2];

i = 1;

%expConditions = [6.02,4,0.5,100] %change this for experimental condition set
ee=1.6022E-19;
me=9.1095E-31; % mass of electron
MO=1836*16*me; % mass of an Oxygen atom
MO2=2*MO; % mass of an Oxygen Molecule
Tg=0.052; % 600K in volts
Ti=Tg;
gammaO2m=0.007;% wall recombination rate of meta-stable Oxygen
R=0.08; % reactor radius
L=0.075; % reactor length
l_p=L/2; % half length
area=2*pi*R*(R+L); % total surface area
volume=pi*R^2*L; % reactor volume
Efactor_O=2+0.5*(1+log(MO/(2*pi*me))); % (E_e+E_i_O)/Te
Efactor_O2=2+0.5*(1+log(MO2/(2*pi*me))); % (E_e+E_i_O2)/Te
% Loading all cross-section data to calculate Kel & Ec
load o2cross.txt -ASCII;
EnergyO2 = o2cross(:,2);
sigmaO2 = o2cross(:,3);
sigO2 = sigmaO2 * 1e-20;
load Ocross.txt -ASCII;
EnergyO = Ocross(:,1);
sigmaO = Ocross(:,2);
sigO = sigmaO * 1e-20;
% power input
Pabs=expConditions(i,4); % total absorbed power in watts [adjustable] 
pabs=Pabs/(ee*volume);
% starting pressures in mTorr (180W)
% atomic oxygen surface recombination rate
gammaO=expConditions(i,3);
p=expConditions(i,1);
Qsccm=expConditions(i,2);
QtorrLit=Qsccm/79.05; % sccm to Torr-Liter/sec
Qmolec=4.483e17*Qsccm; % sccm to molecules/sec
Kpump=2*QtorrLit/(p*volume); % Pumping Rate coefficient
ng0=3.3E19*p*0.026/Tg; % m^-3
ng0_cm=ng0*1e-6; % cm^-3
scat_Xsec=7.5e-19; % elastic scattering cross-section for Oxygen in m^2
% heavy particle reaction rates
Krec=2.6E-14*sqrt(0.026/Tg);
Krec2=2.6E-14*sqrt(0.026/Tg);
Krec3=4.0E-14*sqrt(0.026/Tg);
Krec4=3.3E-17;
Kdet=1.6E-16;
Kch=2.0E-17*sqrt(0.026/Tg);
%
t0=0;
tf=90;
nO2plusbar0=2E16;
nOplusbar0=1E16;
nOminusbar0=3E15;
nObar0=2E17;
nO2mbar0=0.01*ng0;
nO20=ng0-nObar0-nO2mbar0;
Te0=2;
pe0=1.5*(nO2plusbar0+nOplusbar0-nOminusbar0)*Te0;
%v0=[nO20 nObar0 nO2plusbar0 nOplusbar0 nOminusbar0 nO2mbar0 pe0];

% Initial guess to the system
v0=[nO20 nO2plusbar0 nOplusbar0 nOminusbar0 nObar0 nO2mbar0 pe0]/1e19;

%vsoln = fsolve(@steadyState_new_oxys_diss_iz, v0);