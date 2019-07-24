%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-17

function [wa, da] = stark (Ne, Te, L, P, C)
R_Debye=6^(1/3)*C.pi^(1/6)*C.e*(4*C.pi*C.epsilon*C.kB)^(-0.5)*Te^(-0.5)*Ne^(1/6);  %R_Debye:mean distance between ions/Debye radius  
w=spline(L.Tstark, L.w_range, Te)*(Ne/P.Neref);              % w: Electron impact half width in Angstroms 
d=spline(L.Tstark, L.dw_range, Te)*w;                        % d: Electron impact half width in Angstroms
alpha=spline(L.Tstark, L.alpha_range, Te)*(Ne/P.Neref)^0.25; % alpha: Ion broadening factor
wa=(1+1.75*alpha*(1-0.75*R_Debye))*w;                        % Stark half width half maximum in Angstroms
da=(d/w+2*alpha*(1-0.75*R_Debye))*w;                         % Stark shift in Angstroms
%wa=1; da=0;
end