function [Vx Vy] = calcFlux(density,Te,current)
%inputs = density, electron Te, distribution? 
%output = flux
%need to make this universal for any number of ions and neutrals
me=9.1095E-31;
MO=1836*16*me; % mass of an Oxygen atom
MO2=2*MO; % mass of an Oxygen Molecule
ee=1.6022E-19;
density = [10e+15 10e+15 10e+15 10e+15];
nO = density(1);
nO2m = density(2);
nOplus = density(3);
nO2plus = density(4);
final_uB_O2=sqrt(ee*Te/MO2);
final_uB_O=sqrt(ee*Te/MO);
fluxVz_O2 = nO2m*sqrt((8*ee*Te/(pi*MO2)))/4;
fluxVz_O = nO*sqrt((8*ee*Te/(pi*MO)))/4;
fluxVx_O2 = nO2m*sqrt(8*ee*Te*pi/(2*MO2))/2; %cm^2/s
fluxVx_O = nO*sqrt(8*ee*Te*pi/(2*MO))/2; %cm^2/s
nO2plus_flux=nOplus*final_uB_O2*1e-4;
nOplus_flux=nOplus*final_uB_O*1e-4;
fluxVx = [fluxVx_O2 fluxVx_O 0 0 0]
fluxVz = [fluxVz_O2 fluxVx_O nO2plus_flux nOplus_flux]
%calc horizontal velocity
Vx = calcEtchRate(Te,fluxVx,current)
%calc vertical velocity
Vy = calcEtchRate(Te,fluxVz,current)

end