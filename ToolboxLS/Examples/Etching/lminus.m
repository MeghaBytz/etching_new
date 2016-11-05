function lmin=lminus(alpha)
global l_p Rlambda Rrec
Rlam=Rlambda;
Rrc=Rrec;
beta1=max((heff(alpha)-alpha*Rlam)/Rrc./Ffn(alpha), 0);
beta2=max(1-sqrt(alpha*Rlam./heff(alpha)), 0);
lmin=l_p*(1-(beta1.^(-3) + beta2.^(-3)).^(-1/3)); 