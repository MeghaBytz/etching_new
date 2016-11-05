function etchRate = calcEtchRate(Te,fluxes, current)
global noUnknowns
global numberOfIons
global numberOfRadicals
%Surface Kinetics model
%determine free surface sites
%define threshold energy for sputtering and desorption based on Kim et al
%MgO values
epsS = 45; %eV 
epsD = 8; %eV
%default to values of MgO but need to make this generalized
density = 3.56e+3; %kg/m^3
eps = 5.2*Te
M = 40.3;
Na = 6.02e+23;

A = current(1:numberOfIons)
B  =  current(numberOfIons+1:numberOfIons*2);
rxnProb = current(numberOfIons*2+1:numberOfIons*2+numberOfRadicals)
numberOfIons*2+numberOfRadicals+1
length(current)
SCl = current(numberOfIons*2+numberOfRadicals+1:numberOfIons*2+numberOfRadicals*2+1)

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

for i = 1:numberOfRadicals
    chemicalEtchTerm = rxnProb(i)*SCl(i)*radicalFlux(i) + chemicalEtchTerm;
end

surfaceCoverage = chemicalEtchTerm/(desorption+chemicalEtchTerm);

for i = 1:numberOfRadicals
    chemicalEtch = (1-surfaceCoverage)*rxnProb(i)*radicalFlux(i) + chemicalEtch;
end

for i = 1:numberOfRadicals
    sputtering = (1-surfaceCoverage)*sputteringYield(i)*ionFlux(i) + sputtering;
end
%EtchRate = (1-RxnProb*FluxCl/(RxnProb*FluxCl + Yd*TotalPosFlux))*(RxnProb*FluxCl+Ys*FluxAr_pos);
etchRate = chemicalEtch + sputtering;
etchRate = etchRate*6e+11*M/(Na*density)+current(end)*randn; %nm/min
%add unit conversions and velocity terms to make -- fluxes in units of 
end


    


