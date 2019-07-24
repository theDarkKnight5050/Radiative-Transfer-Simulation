%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Computes ion and neutral number density along with volumetric electron density assuming LTE

function [Ne, N1, N2, Zci] = saha(beta, No, L, C)

[Zci Zcn]=partition(L,(C.kB*beta)^-1);

%solve for Ne taking into account ionization depression
const=(2*Zci/Zcn)*((2*C.pi*C.me/beta/(C.h^2))^1.5)*exp(-L.Ei*beta);
Ne=(-const+sqrt(const^2+4*No*const))/2; %Saha equation for electron density (quadratic)
if Ne>1e20 %only starts taking into account depression at 10^20m^(-3)
    while true 
        E_id=spline(L.Nedep,L.dE1,log10(Ne));
        const=(2*Zci/Zcn)*(2*C.pi*C.me/beta/C.h^2)^1.5*exp(-(L.Ei-E_id*C.e)*beta);
        Ne=(-const+sqrt(const^2+4*No*const))/2;
        if abs(E_id-spline(L.Nedep,L.dE1,log10(Ne)))<1e-2, break; %break when improvement becomes negligible
        else, E_id=spline(L.Nedep,L.dE1,log10(Ne));
        end 
        %coeff2=(Z_On^2/Z_O2)*(1/2)^1.5*(2*pi*m_O*k_B*T_e/h^2)^1.5*exp(-E_d/k_B/T_e);        
        %fun=@(x) [x(3)^2/x(2)-coeff1; x(2)^2/x(1)-coeff2; 2*x(1)+x(2)+x(3)-2.67e25*0.21*2];
        %x0=[1e16; 1e16; 1e16];
        %[N_O2, N_O, N_e]=fsolve(fun,x0);
    end
end

N2=(No-Ne)*(L.g2a/Zcn)*exp(-L.E2*beta);
N1=N2*(L.g1a/L.g2a)*exp(L.E21*beta); %Boltzam distribution

end
