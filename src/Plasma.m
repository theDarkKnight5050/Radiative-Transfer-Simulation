%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Class storing plasma constants
classdef Plasma
  properties
      Rp=0.00015;  %radius of the plasma (for parabola- hard cut-off, for Gaussian - core)
      Ro=0.00025;  %outer boundary when using a Gaussian temperature profile
      Temax;       %peak plasma temperature
      Temin=2000;  %minimum plasma temperature
      No=2.67e25;  %atom density (m^-3)- we assume plasma has not had time to expand
      Neref=1e22;  %reference density from Griem (m^-3);
      Wp;          %plasma frequency=Wp*ne^0.5
  end
  methods
    function plasma=Plasma(C)
      plasma.Temax=11604*input('Input max temperature in electron volts (minimum of 0.172): ');
      plasma.Wp=(C.e^2/C.me/C.epsilon)^0.5;
    end
  end
end