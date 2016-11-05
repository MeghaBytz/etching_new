function vmin=vol_minus(alpha)
vmin=4*pi/3.*(rminus(alpha).^2 - rminus(alpha).*lminus(alpha)*2/3 + (lminus(alpha).^2)/6).*lminus(alpha);