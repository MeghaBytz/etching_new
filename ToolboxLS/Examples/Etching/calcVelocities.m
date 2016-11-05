function [Vx Vy] = calcVelocities(density,Ti,Te,current)
%inputs = density, electron Te, distribution? 
%output = flux
%need to make this universal for any number of ions and neutrals
me=9.1095E-31;
MO=1836*16*me; % mass of an Oxygen atom
MO2=2*MO; % mass of an Oxygen Molecule
ee=1.6022E-19;
nO = density(1);
nO2m = density(2);
nOplus = density(3);
nO2plus = density(4);
final_uB_O2=sqrt(ee*Te/MO2);
final_uB_O=sqrt(ee*Te/MO);
fluxVz_O2 = nO2m*sqrt((8*ee*Ti/(pi*MO2)))/4;
fluxVz_O = nO*sqrt((8*ee*Ti/(pi*MO)))/4;
fluxVx_O2 = nO2m*sqrt((8*ee*Ti/(pi*MO2)))/4; %cm^2/s %check Maxwellian assumption for ions
fluxVx_O = nO*sqrt((8*ee*Ti/(pi*MO)))/4;
fluxVx_Oplus = nOplus*final_uB_O2/(2*pi)*1e-4;
fluxVx_O2plus = nOplus*final_uB_O2/(2*pi)*1e-4;%cm^2/s
nO2plus_flux=nOplus*final_uB_O2*1e-4;
nOplus_flux=nOplus*final_uB_O*1e-4;
fluxVx = [ fluxVx_O2plus   fluxVx_Oplus  fluxVx_O2 fluxVx_O]
assignin('base', 'fluxVx', fluxVx)
assignin('base', 'Te', Te)
fluxVz = [ nO2plus_flux nOplus_flux fluxVz_O2 fluxVz_O]
assignin('base', 'fluxVz', fluxVz)
%calc horizontal velocity
Vx = calcEtchRate(Te,fluxVx,current)/1000;
%calc vertical velocity
Vy = calcEtchRate(Te,fluxVz,current)/1000; %divide by 1000 for um and divide by 10 to normalize for grid
if Vx<0
    Vx = 0;
end
if Vy<0
    Vy = 0;
end
assignin('base', 'Vy', Vy)
    assignin('base', 'Vx', Vx)
end