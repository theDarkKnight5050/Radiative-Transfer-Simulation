%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Calculates radiative transfer coefficients

function [eps, kap] = bremss (n1, n2, L, C, ne, Te)

Pa=(2*C.h*C.co^2/((L.lama*1e-10)^5))/(exp(L.E21/C.kB/Te)-1); %Planck body source function 
eps=L.E21/(4*C.pi)*n2*L.A21a;
%kap=L.E21/(4*C.pi)*(C.co/L.lama^2)*(n1*L.B12a-n2*L.B21a);
%kap=0;
kap=eps/Pa;
%kap=1.37E-27*ne^2*L.lama^3/Te^0.5*(1-exp(-L.E21/C.kB/Te));

endfunction
