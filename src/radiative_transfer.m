%This program computes the analytical solution for the 1D radiation transfer equation for a spectral
%line from a cylindrically symmetric source - M.A. Cappelli June 2019
clear all;
close all;
disp("This program computes the analytical solution for the 1D radiation transfer equation for a spectral line from a cylindrically symmetric source");
disp("Author: M.A. Cappelli");

%fundamental constants
C=Constants;

% Defines parameters of the current line
disp("Computing line properties...");
L=Line(C);

%% Deals with the geometry and plasma conditions
disp("Loading plasma conditions...");
P=Plasma;
prompt="Input the y-position of interest (must be less than %d if parabolic Te, less than %d if Gaussian): ";
ypos=input(sprintf(prompt, P.Rp, P.Ro)); 
%Gaussian Temp Profile
xmin=-(P.Rp^2-ypos^2)^.5;
xmin=-(P.Rp^2-ypos^2)^.5;
xmax=(P.Rp^2-ypos^2)^.5;
%Parabolic Temp Profie
%xmin=-(P.Ro^2-ypos^2)^.5;
%xmax=(P.Ro^2-ypos^2)^.5;
xx=xmin:(xmax-xmin)/40:xmax;

numPoints=input(prompt="Input the desired number of data points: ");
x=1000/numPoints;
np=-500:x:500;  %advancement of wavelength in Angstroms
intgrl=zeros(1,length(np));
disp("Processing each spectral element (each . represents 5 dlambdas)");
count=0;
colors=varycolor(numPoints+1);
for i=1:length(np) %spectral frequency interval
    count++;
    if(count==5)
        printf(".");
        count=0;
    endif
    
    npx=np(i);
    epsval=[];
    kapval=[];
    epsL=[];
    kapL=[];
    
    for k=1:length(xx)
        x=xx(k);
        Rxy=(x^2+ypos^2)^.5;
        
        %%Temperature calculation
        %Parabolic Temp
        %Te=Temin+Temax*(1-0.5*Rxy^2/Rp^2); 
        %Gaussian Temp
        Te=(P.Temax-P.Temin)*exp(-(Rxy^2)/(P.Rp^2))+P.Temin;
        beta=1/(C.kB*Te);
        Tex(k)=Te;

        %% Computes the LTE model
        % Saha Equation
        [ne n1 n2 Zci]=saha(beta, P.No, L, C);
        Nex(k)=ne;
        N2x(k)=n2;
        N1x(k)=n1;
        % Emission and Absorption        
        [eps kap]=bremss(n1, n2, Zci, npx, L, P, C, ne, Te);       
        epsL(i,k)=eps;
        kapL(i,k)=kap;
        
        %pcoeff=8;  %coefficient of polynomial fit.
        % plasma geometry is a cylinder
        % %Fitting functions to the volume emission and absorption accross slabs
        % epsfit=polyfit(xx, epsL, pcoeff);
        % kapfit=polyfit(xx, kapL, pcoeff);
    end
    figure(1)
    plot(xx,(epsL(i,:)), 'Color',colors(i,:));
    % %plot(xx,epsL,'r-', xx, epsval(i,:), 'ro');
    ylabel('Volume Emission Coefficient (log)');
    xlabel('position');
    hold on
    % 
    figure(2)
    plot(xx,(kapL(i,:)), 'Color',colors(i,:));
    %plot (xx,kapL,'r-', xx, kapval(i,:), 'ro');
    ylabel('Volume Absorption Coefficient (log)');
    xlabel('position');
    hold on

    for k = 1:length(xx) 
        for kk=k:length(xx)
            index=kk+1-k;
            xxk(index)=xx(kk);
            intgndkap(index)=kapL(i,kk); 
        end
        for kk=k:length(xx)
            intglkap(kk)=trapz(xxk,intgndkap);
        end
        integrand(k)=epsL(i,k)*exp(-intglkap(k));         
    end
    intgrl(i)=trapz(xx,integrand);
end
printf("\n"); 

input("Press [enter] to graph results");
disp("Graphing results...");

figure(3)
semilogy(xx,N1x,'r-');
ylabel('Lower State Density');
xlabel('position');
hold on

figure(4)
semilogy(xx,N2x,'r-');
ylabel('Upper State Density');
xlabel('position');
hold on

figure(5)
plot(xx,Tex,'r-');
ylabel('Temperature');
xlabel('position');
hold on

figure(6)
semilogy(xx,Nex,'r-');
ylabel('Electron Density');
xlabel('position');
hold on

figure(7)
plot(np,intgrl,'r-');
ylabel('spectral radiance');
xlabel('wavelength(nm)'); 