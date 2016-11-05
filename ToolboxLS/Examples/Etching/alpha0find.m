function G=alpha0find(alpha)
global volume alphabar
alpha_bar=alphabar;
volm=volume;
G=alpha.*vol_minus(alpha) - alpha_bar*volm; 