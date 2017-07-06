function [ pps,TTs,p0p0s,VVs,rhorhos,fLD ] = FannoM( M )
%Input:M, Outputs: p/pstar, T/Tstar, p0/p0star, V/Vstar, rho/rhostar,4fLstar/D
%   Written: Chris Kniffin, Spring 2016
g = 1.4;

fLD = (1-M^2)/(g*M^2)+(g+1)/(2*g)*log(M^2/((2/(g+1)*(1+(g-1)/2*M^2))));
pps = 1/M*1/sqrt((2/(g+1)*(1+(g-1)/2*M^2)));
rhorhos = 1/M*sqrt(2/(g+1)*(1+(g-1)/2*M^2));
TTs = 1/(2/(g+1)*(1+(g-1)/2*M^2));
VVs = M/sqrt(2/(g+1)*(1+(g-1)/2*M^2));
p0p0s = 1/M*(2/(g+1)*(1+(g-1)/2*M^2))^((g+1)/(2*(g-1)));

end

