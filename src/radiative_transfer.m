%This program computes the analytical solution for the 1D radiation transfer equation for a spectral line from a 
%cylindrically symmetric source - M.A. Cappelli June 2019
clear all;
close all;
disp('This program computes the analytical solution for the 1D radiation transfer equation for a spectral line from a cylindrically symmetric source');
disp('Author: M.A. Cappelli');

% Fundamental constants
C=Constants;

% Defines parameters of the current line
disp('Computing line properties...');
L=Line(C);

%% Deals with the geometry and plasma conditions
disp('Loading plasma conditions...');
P=Plasma(C);
prompt='Input the y-position of interest (must be less than %d if parabolic Te, less than %d if Gaussian): ';
ypos=input(sprintf(prompt, P.Rp, P.Ro)); 
% Gaussian Temp Profile
xmin=-(P.Ro^2-ypos^2)^.5;
xmax=(P.Ro^2-ypos^2)^.5;
% Parabolic Temp Profie
%xmin=-(P.Rp^2-ypos^2)^.5;
%xmax=(P.Rp^2-ypos^2)^.5;
xx=xmin:(xmax-xmin)/40:xmax;

numPoints=input('Input the desired number of data points: ');
width=input('Input the desired window size (in Angstroms): ');
x=2*width/numPoints;
np=-width:x:width;  %advancement of wavelength in Angstroms
intgrl=zeros(1,length(np));
disp('Processing each spectral element (each . represents 5 dlambdas)');
count=0;
colors=varycolor(numPoints+1);
for i=1:length(np) %spectral frequency interval
    count=count+1;
    if(count==5)
        fprintf('.');
        count=0;
    end
    
    npx=np(i);
    epsL=[];
    kapL=[];
    ratioL=[];
    
    for k=1:length(xx)
        x=xx(k);
        Rxy=(x^2+ypos^2)^.5;
        
        %% Temperature calculation
        % Parabolic Temp
        %Te=P.Temin+(P.Temax-P.Temin)*(1-0.5*Rxy^2/P.Rp^2); 
        % Gaussian Temp
        Te=(P.Temax-P.Temin)*exp(-(Rxy^2)/(P.Rp^2))+P.Temin;
        beta=1/(C.kB*Te);
        Tex(k)=Te;

        %% Computes the LTE model
        % Saha Equation
        [ne n1 n2 Zci]=saha(beta, P.No, L, C);
        Nex(k)=ne;
        N1x(k)=n1;
        N2x(k)=n2;
        % Emission and Absorption        
        [eps kap ratio]=bremss(n1, n2, Zci, npx, L, P, C, ne, Te);       
        epsL(i,k)=eps;
        kapL(i,k)=kap;
        ratioL(i,k)=ratio;       
    end
    
    figure(1)
    plot(xx,log10(epsL(i,:)),'Color',colors(i,:));
    ylabel('Volume Emission Coefficient (log)');
    xlabel('position');
    hold on
    % 
    figure(2)
    plot(xx,log10(kapL(i,:)),'Color',colors(i,:));
    ylabel('Volume Absorption Coefficient (log)');
    xlabel('position');
    hold on
    %
    figure(3)
    plot(xx,log10(ratioL(i,:)),'Color',colors(i,:));
    ylabel('Line-Continuum Emission Ratio');
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
fprintf('\n'); 

input('Press [enter] to graph results');
disp('Graphing results...');

figure(4)
plot(xx,Tex,'r-');
ylabel('Temperature');
xlabel('position');
hold on

figure(5)
semilogy(xx,Nex,'r-');
ylabel('Electron Density');
xlabel('position');
hold on

figure(6)
semilogy(xx,N1x,'r-');
ylabel('Lower State Density');
xlabel('position');
hold on

figure(7)
semilogy(xx,N2x,'r-');
ylabel('Upper State Density');
xlabel('position');
hold on

figure(8)
plot(np,(intgrl),'r-');
ylabel('spectral radiance');
xlabel('wavelength(A)'); 