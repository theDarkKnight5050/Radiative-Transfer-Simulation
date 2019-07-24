%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Class storing line constants 
%% Current line: O7771.94
classdef Line
    properties
        A21a=3.69e7;  %einstein spontaneous emission coefficient (s^-1)
        lama=7771.94; %line a center wavelength (note Oxygen is a triplet) in Angstroms
        g1a=5;        %lower state generacy
        g2a=15;       %upper e degeneracy - note that for O we are lumping the three states together
        Ei;           %ionization energy (J) of the oxygen atom (not reduced by shielding)
        E2;           %upper state energy  (J)
        E21;          %transition energy (J)
        E1;           %lower state energy (J)
        B21a;         %Einstein Stimulated Emission (See Wikipedia - Einstein Coefficients) given in wavelength units for spectral radiance [B]=M^(-1)LT^2
        B12a;         %Einstein Stimulated Absorption
        
        %Constants related to the partition function as determined from M. Capitelli et al.Tables of Internal Partition 
        %Functions and Thermodynamic Properties of High-Temperature Mars-Atmosphere Species
        Tpart; %Associated temperature
        Yn;    %Neutral species
        Y1;    %First ionization
        Y2;    %Second ionization
        
        Nedep; %Associated electron density
        dE1;   %Ionization depression
    
        %Recalling electron impact half width values (Griem)
        Tstark=[2500 5000 10000 20000 40000 80000];                %Associated temperature
        w_range=[1.99e-2 2.48e-2 3.27e-2 4.43e-2 5.66e-2 6.49e-2]; %electron impact half width
        dw_range=[1.768 1.676 1.421 1.076 0.783 0.585];            %shift over width ratio
        alpha_range=[0.041 0.035 0.028 0.023 0.019 0.017];         %ion broadening parameter
    end
    methods
        function line=Line(C)
            line.Ei=(13.62)*C.e; 
            line.E2=86628.7*C.h*100*C.co; %https://physics.nist.gov/PhysRefData/Handbook/Tables/oxygentable5.htm
            line.E21=(C.h*C.co/(line.lama*1e-10)); 
            line.E1=line.E2-line.E21; 
            line.B21a=line.A21a*((line.lama*1e-10)^5/(2*C.h*C.co^2));  
            line.B12a=line.g2a*line.B21a/line.g1a; 
            
            M=dlmread('../data/partition_function.csv',',',1,0);
            line.Tpart=M(:,1); 
            line.Yn=M(:,2);    
            line.Y1=M(:,3);    
            line.Y2=M(:,4);    
            
            M=dlmread('../data/ionization_depression.csv',',',1,0);
            line.Nedep=log10(M(:,1)); 
            line.dE1=M(:,2);          
        end
   end
end