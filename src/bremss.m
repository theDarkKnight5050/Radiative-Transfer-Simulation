%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Calculates radiative transfer coefficients

function [eps, kap, ratio] = bremss (n1, n2, Zi, pos, L, P, C, ne, Te)

eps=0;
kap=0;

lama=(L.lama+pos)*1E-10;
[wa da]=stark(ne, Te, L, P, C);
E=C.h*C.co/lama;

Linshape=(1/C.pi/wa)/(((pos-da)/wa)^2+1)*10^10;
Pa=2*C.h*C.co^2/lama^5/(exp(E/C.kB/Te)-1); %Planck body source function

ratio=L.E21/(4*C.pi)*n2*L.A21a*Linshape;
eps=eps+ratio;                                                              %bound-bound emission
eps=eps+1.63E-43*ne^2*lama^(-2)*Te^(-0.5)*L.g2a/Zi*(1-exp(-E/C.kB/Te))*100; %free-bound emisson
eps=eps+1.63E-43*ne^2*lama^(-2)*Te^(-0.5)*exp(-E/C.kB/Te)*100;              %free-free emission
if C.co/lama>P.Wp*ne^0.5
  kap=kap+eps/Pa;
end
ratio=ratio/(eps-ratio);
%kap+=L.E21/(4*C.pi)*(n1*L.B12a-n2*L.B21a)*Linshape;
%kap+=1.37E-27*ne^2*lama^3*Te^(-0.5)*(1-exp(-L.E21/C.kB/Te));
%kap+=0;
%eps+=kap*Pa;

end