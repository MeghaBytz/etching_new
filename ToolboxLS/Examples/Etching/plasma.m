function [Vx Vy] = plasma(current,i)
global ee me MO MO2 ng Tg Ti l_p gammaO gammaO2m
global Efactor_O2 Efactor_O EnergyO2 sigO2 EnergyO sigO
global Krec Krec2 Krec3 Krec4 Kdet Kch Rlambda hl0 Rrec alphabar
global pabs R area volume QtorrLit Qmolec Kpump scat_Xsec
global numberOfIons
global numberOfRadicals
global expConditions
global viewFactor
global noUnknowns

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
length(expConditions)
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
v0=[nO20 nO2plusbar0 nOplusbar0 nOminusbar0 nObar0 nO2mbar0 pe0];
[t v]=ode23s('new_oxys_diss_iz', [t0 tf], v0);
nO2=v(:,1);
nO2plusbar=v(:,2);
nOplusbar=v(:,3);
nOminusbar=v(:,4);
ne0=nO2plusbar+nOplusbar-nOminusbar;
nObar=v(:,5);
nO2mbar=v(:,6);
Te=v(:,7)./(1.5.*ne0);
% getting the data at the final equilibrium
final_nO2=nO2(end);
n_O2=final_nO2*1e-6; % in cm^-3
final_nO2plusbar=nO2plusbar(end);
n_O2plus_bar=final_nO2plusbar*1e-6; % in cm^-3
final_nOplusbar=nOplusbar(end);
n_Oplus_bar=final_nOplusbar*1e-6; % in cm^-3
final_nOminusbar=nOminusbar(end);
n_Ominus_bar=final_nOminusbar*1e-6; % in cm^-3
final_ne0=ne0(end);
n_e0=final_ne0*1e-6; % in cm^-3
n_e0_2=n_e0/1e10; 
final_nObar=nObar(end);
n_O_bar=final_nObar*1e-6; % in cm^-3
final_nO2mbar=nO2mbar(end);
n_O2m_bar=final_nO2mbar*1e-6; % in cm^-3
final_ng=final_nO2+final_nObar+final_nO2mbar;
n_g=final_ng*1e-6; % in cm^-3
final_Te=Te(end);
final_Ti=Ti;
final_Ec_O2=o2_ec(final_Te);
final_Ec_O=o_ec(final_Te);
final_p=final_ng/3.3e19*(Tg/0.026); % final pressure in mTorr
% ratio for density weighting
nOp_ratio=final_nOplusbar/(final_nO2plusbar+final_nOplusbar);
% plotting the results
% figure(ii)
% subplot(7,1,1)
% plot(t,nO2plusbar,t, ne0,'--')
% ylabel('n_{{O_2}^+}, n_{e0} (m^{-3})')
% axis([0 inf 0 final_nO2plusbar*1.5])
% title(['Flowrate=',num2str(Qsccm),'sccm, P_{abs}=',num2str(Pabs),...
%  'W, p_0=',num2str(p),'mTorr, p_f=',num2str(round(final_p)),...
%  'mTorr'])
% subplot(7,1,2)
% plot(t,nOplusbar)
% ylabel('n_{O^+} (m^{-3})')
% axis([0 inf 0 final_nOplusbar*1.5])
% subplot(7,1,3)
% plot(t,nOminusbar)
% ylabel('n_{O^-} (m^{-3})')
% axis([0 inf 0 final_nOminusbar*1.5])
% subplot(7,1,4)
% plot(t,nObar)
% ylabel('n_O (m^{-3})')
% axis([0 inf 0 final_nObar*1.5])
% subplot(7,1,5)
% plot(t,nO2mbar)
% ylabel('n_{{O_2}*} (m^{-3})')
% axis([0 inf 0 final_nO2mbar*1.5])
% subplot(7,1,6)
% plot(t,nO2)
% ylabel('n_{O_2} (m^{-3})')
% axis([0 inf 0 final_nO2*1.5])
% subplot(7,1,7)
% plot(t,Te)
% xlabel('t (sec)')
% ylabel('T_e (m^{-3})')
% axis([0 inf 0 final_Te*1.5])
% caculating values for O neutral & O2m
lambda=1/(final_ng*scat_Xsec); % lambda in m
vbarO=sqrt(8*ee*Ti/(pi*MO)); % average thermal velocity of O neutral
vbarO2m=sqrt(8*ee*Ti/(pi*MO2)); % average thermal velocity of O2m
DO=ee*Tg*lambda/vbarO/MO; % Diffusion coefficient of O neutral
DO2m=ee*Tg*lambda/vbarO2m/MO2; % Diffusion coefficient of O2m
dO=sqrt(4*DO*l_p*(2-gammaO)/vbarO/gammaO + l_p^2);
dO2m=sqrt(4*DO2m*l_p*(2-gammaO2m)/vbarO2m/gammaO2m + l_p^2);
hAO=1/(1 + l_p*vbarO*gammaO/4/DO/(2-gammaO));
hAO2m=1/(1 + l_p*vbarO2m*gammaO2m/4/DO2m/(2-gammaO2m)) ;
vol_O=volume*(1 - l_p^2/(3*dO^2))*(1 - (2/3)*l_p^3/(R*dO^2)...
 + l_p^4/(6*R^2*dO^2));
vol_O2m=volume*(1 - l_p^2/(3*dO2m^2))*(1 - (2/3)*l_p^3/(R*dO2m^2)...
 + l_p^4/(6*R^2*dO2m^2));
vr_O=vol_O/volume;
vr_O2m=vol_O2m/volume;
% calculation of the final alpha0
Tplusf=Ti;
Tminusf=Ti;
gamma_plusf=final_Te/Tplusf;
gamma_minusf=final_Te/Tminusf;
etaf=2*Tplusf/(Tplusf+Tminusf);
hl0=0.86/sqrt(3+etaf*l_p/lambda);
Rlambda=sqrt(2*pi/gamma_plusf)*lambda/l_p/etaf;
final_uB_O2=sqrt(ee*final_Te/MO2);
final_uB_O=sqrt(ee*final_Te/MO);
% density-weighted Bohm velocity in m/s
final_uB_dw=final_uB_O2*(1-nOp_ratio)+final_uB_O*nOp_ratio;
% density-weighted Recombination rate coefficient
Krec_dw_f=(Krec+Krec2)*(1-nOp_ratio)+Krec3*nOp_ratio;
Rrec=Krec_dw_f*final_ne0*l_p/final_uB_dw;
alphabar=final_nOminusbar/final_ne0;
alpha0=fzero('alpha0find',[1e-2 1e2]); % final alpha0
volumeminus=vol_minus(alpha0);
final_nOminus=final_nOminusbar*volume/volumeminus;
final_nplus=final_nOminus + final_ne0;
final_nOplus=final_nplus*nOp_ratio;
final_nO2plus=final_nplus-final_nOplus;
final_nO=final_nObar*volume/vol_O;
final_nO2m=final_nO2mbar*volume/vol_O2m;
n_O2plus=final_nO2plus*1e-6; % in cm^-3
n_Oplus=final_nOplus*1e-6; % in cm^-3
n_Ominus=final_nOminus*1e-6; % in cm^-3
n_O=final_nO*1e-6; % in cm^-3
n_O2m=final_nO2m*1e-6; % in cm^-3
% calculation of final hl factor
% density weighted positive ion mass
Mplus_dw=MO2*(1-nOp_ratio)+MO*nOp_ratio;
nstarf=15/56*sqrt(8*ee*Tplusf/pi/Mplus_dw)*(etaf^2)/(Krec_dw_f*lambda);
hpar2f=hl0*1/(1+alpha0);
hpar1f=1/(gamma_minusf^0.5 + (gamma_plusf^0.5)...
 *(etaf*l_p/sqrt(2*pi)/lambda))*(alpha0/(1+alpha0));
hflat1f=1/(gamma_minusf^0.5 + (gamma_plusf^0.5)*(nstarf^0.5)...
 /(final_nOminus^0.5));
%overall hl factor by "linear ansatz"
final_hl=hpar2f + hpar1f + hflat1f;
% size of EN core (lminus & rminus)
final_lminus=lminus(alpha0); %in m
Lminus=final_lminus*1e2; % in cm
final_rminus=R - l_p + final_lminus; %in m
Rminus=final_rminus*1e2; % in cm
lminus_over_l_p=final_lminus/l_p;
vratio=volumeminus/volume;
vol_rec=(2*pi*final_lminus/(1+alpha0))*(8/15*alpha0*...
 (final_rminus^2 - 14/15*final_rminus*final_lminus...
 + 4/15*final_lminus^2) + 2/3*(final_rminus^2 - 2/3*final_rminus*final_lminus + 1/6*final_lminus^2)); 
vrec_ratio=vol_rec/volume;
% ion flux
% final_nplus=final_nOminus + final_ne0;
% final_nOplus=final_nplus*nOp_ratio;
% final_nO2plus=final_nplus-final_nOplus;
% final_nO=final_nObar*volume/vol_O;
% final_nO2m=final_nO2mbar*volume/vol_O2m;
dissociation_rate = final_nObar/final_ng;

%need to make these calculations general
nplus_flux=final_hl*final_nplus*final_uB_dw*1e-4 % in /cm^2/s
nOplus_flux=final_hl*final_nOplus*final_uB_dw*1e-4 % in /cm^2/s
nO2plus_flux=final_hl*final_nplus*final_uB_dw*1e-4 % in /cm^2/s
nO_flux=hAO*final_nO*vbarO/4*1e-4 % in /cm^2/s
nO2m_flux=hAO*final_nO2m*vbarO2m/4*1e-4; % in /cm^2/s
flux_ratio=nO_flux/nplus_flux;
lambda_cm=lambda/1e-2; % in cm

densities = [ final_nO2m final_nO final_nOplus final_nO2plus ]
assignin('base', 'densities', densities)
[Vx Vy] = calcVelocities(densities, final_Ti,final_Te,current); %divide by 100 to normalize for grid that goes crom -100 to 100
%fluxes = [nO2plus_flux nOplus_flux nO_flux nO2m_flux]

%Calculate etch rate
%etchRate = calcEtchRate(final_Te,fluxes, current)
% save results
% allresults(:,ii)=[p;final_p;Pabs;ng0_cm;n_g;n_O2plus;n_Oplus;...
%  n_Ominus;n_e0;n_e0_2;n_O;n_O2m;n_O2plus_bar;n_Oplus_bar;...
%  n_Ominus_bar;n_O_bar;n_O2m_bar;n_O2;final_Te;final_Ec_O2;...
%  final_Ec_O; alpha0;alphabar;gamma_plusf;final_hl;...
%  lminus_over_l_p;Rminus;Lminus;vratio;vrec_ratio;...
%  lambda_cm;dissociation_rate;nplus_flux;nO_flux;flux_ratio];

% save results to a file
% filename=[num2str(Pabs),'W_results.txt']
% save(filename,'allresults','-ASCII','-double'); 