function [ M,pps ] = FannoFLDsub( fLD )
%Input: 4fL/D Outputs: M (subsonic), p/pstar
%   Written: Chris Kniffin, Spring 2016

Ml = 0;
Mu = 1;
M = .1;
k = 0;
kmax = 1000;
g = 1.4;
fLDc = 0;

while abs(fLD-fLDc) > eps && k < kmax
    fLDc = (1-M^2)/(g*M^2)+(g+1)/(2*g)*log(M^2/((2/(g+1)*(1+(g-1)/2*M^2))));
    if fLDc > fLD
        Ml = M;
        M  = (Mu+M)/2;
    else
        Mu = M;
        M  = (Ml+M)/2;
    end
    k = k+1;
end

pps = 1/M*1/sqrt((2/(g+1)*(1+(g-1)/2*M^2)));