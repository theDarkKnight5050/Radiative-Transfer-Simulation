function linshape = voigt (wl, wg, dl, a, pos)

wv=wl/2+((wl/2)^2+wg^2)^0.5;
linshape=(1-wl/wv)exp(-4*ln(2)*((dl-pos)/wv)^2)+wl/wv/(1+4*((dl-pos)/wv)^2);

end
