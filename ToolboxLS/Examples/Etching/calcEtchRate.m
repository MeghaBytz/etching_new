function etchRate = calcEtchRate(Te,fluxes, current)
global noUnknowns
global numberOfIons
global numberOfRadicals
%Surface Kinetics model
%determine free surface sites
%define threshold energy for sputtering and desorption based on Kim et al
%MgO values
current
epsS = current(11); %eV 
epsD = current(12); %eV
%default to values of MgO but need to make this generalized
density = 1180; %kg/m^3 make density of resist materials
eps = 5.2*Te;
M = 100.12;
Na = 6.02e+23;

A = current(1:2)
B  =  current(3:4);
rxnProb = current(5:6)
SCl = current(7:8);

sputteringYield = A.*(sqrt(eps) - sqrt(epsS));
ionStimulatedDesorption = B.*(sqrt(eps) - sqrt(epsD))

%set negative entries equal to 0
sputteringYield(sputteringYield<0)=0;
ionStimulatedDesorption(ionStimulatedDesorption<0)=0;

%neglecting O2m right now and thermal desorption
%fluxes = [nO2plus_flux nOplus_flux nO_flux nO2m_flux]
ionFlux = fluxes(1:numberOfIons);
radicalFlux = fluxes(numberOfIons+1:length(fluxes));
chemicalEtchTerm = 0; %gamma*radicalFlux
desorption = 0;
chemicalEtch = 0;
sputtering = 0;

for i = 1: numberOfIons
    desorption = ionStimulatedDesorption(i)*ionFlux(i) + desorption;
end
assignin('base', 'radicalFlux', radicalFlux)
assignin('base', 'ionFlux', ionFlux)
assignin('base', 'ionStimulatedDesorption', ionStimulatedDesorption)
assignin('base', 'desorption', desorption)
for i = 1:numberOfRadicals
    chemicalEtchTerm = rxnProb(i)*SCl(i)*radicalFlux(i) + chemicalEtchTerm;
end
assignin('base', 'chemicalEtchTerm', chemicalEtchTerm)
surfaceCoverage = chemicalEtchTerm/(desorption+chemicalEtchTerm);
assignin('base', 'surfaceCoverage', surfaceCoverage)
for i = 1:numberOfRadicals
    chemicalEtch = (1-surfaceCoverage)*rxnProb(i)*radicalFlux(i) + chemicalEtch;
end
assignin('base', 'chemicalEtch', chemicalEtch)
for i = 1:numberOfRadicals
    sputtering = (1-surfaceCoverage)*sputteringYield(i)*ionFlux(i) + sputtering;
end
assignin('base', 'sputtering', sputtering)
%EtchRate = (1-RxnProb*FluxCl/(RxnProb*FluxCl + Yd*TotalPosFlux))*(RxnProb*FluxCl+Ys*FluxAr_pos);
etchRate = chemicalEtch + sputtering;
etchRate = etchRate*6e+11*M/(Na*density);%current(end)*randn; %nm/min
%add unit conversions and velocity terms to make -- fluxes in units of 
end


    


