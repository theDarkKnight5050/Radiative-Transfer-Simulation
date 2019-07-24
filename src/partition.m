%% Author: aiabd <aiabd@LAPTOP-R5RTKBLK>
%% Created: 2019-07-11
%% Compute ion and neutral partition function at some T given polynomial fits

function [Zcn, Zci] = partition (L, Te)
 
Zcn=spline(L.Tpart,L.Yn,Te); %neutral partition function
Zci=spline(L.Tpart,L.Y1,Te); %ion partition function
Zc2=spline(L.Tpart,L.Y2,Te);

end
