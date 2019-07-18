%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Class storing plasma constants
classdef Plasma
  properties(Constant=true)
      Rp = 0.00015;  %radius of the plasma (for parabola- hard cut-off, for Gaussian - core)
      Ro=0.00025; %outer boundary when using a Gaussian temperature profile
      Temax=25000;  %peak plasma temperature
      Temin =2000; %minimum plasma temperature
      No = 2.67e25; %atom density (m^-3)- we assume plasma has not had time to expand
      Neref=1e22; %reference density from Griem (m^-3);
  end
endclassdef
