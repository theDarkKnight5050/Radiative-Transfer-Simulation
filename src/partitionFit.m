%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Fits partition function 7th degree polynomial

function [fi, fn] = partitionFit(Yi, Yn)

warning('off', 'Octave:nearly-singular-matrix');
%X=[1000 5000 10000 15000 20000 25000 30000 36000];
fi=polyfit(X, Yi, 7);  %describes the polynomial for the ionized partition function
fn=polyfit(X, Yn, 7);  %describes the polynomial for the neutral partition function

end
