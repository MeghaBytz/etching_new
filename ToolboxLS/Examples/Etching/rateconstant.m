function K=rateconstant( Energy, sig, Tel );
global ee me
nefnull = (2 * ee .* Energy /me).^(1/2) .* sig;
q0 = Energy.^0.5.*2/sqrt(pi).*(1/Tel)^(3/2).*exp(-Energy/Tel).*nefnull;
qq=0;qqq=0;
for ii=1:1:max(size(Energy))-1
 qq = (q0(ii) + q0(ii+1)) .* (Energy(ii+1) - Energy(ii))/2;
 qqq = qq + qqq;
end
K = qqq; 