%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Calculates radiative transfer coefficients

function [eps, kap] = bremss (n1, n2, Zi, lama, L, C, ne, Te)

eps=kap=0;
Pa=2*C.h*C.co^2/lama^5/(exp(L.E21/C.kB/Te)-1); %Planck body source function 
eps+=L.E21/(4*C.pi)*n2*L.A21a;                                           %bound-bound emission
eps+=1.63E-43*ne^2*lama^(-2)*Te^(-0.5)*L.g2a/Zi*(1-exp(-L.E21/C.kB/Te)); %free-bound emisson
eps+=1.63E-43*ne^2*lama^(-2)*Te^(-0.5)*exp(-L.E21/C.kB/Te);              %free-free emission
kap+=eps/Pa;
%kap+=L.E21/(4*C.pi)*(n1*L.B12a-n2*L.B21a);
%kap+=1.37E-27*ne^2*lama^3*Te^(-0.5)*(1-exp(-L.E21/C.kB/Te));
%kap+=0;
%eps+=kap*Pa;

endfunction
