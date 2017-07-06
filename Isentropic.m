function [ p0p,A_As,T0T,rho0rho ] = Isentropic( M )
%Outputs: p0/p, A/As, T0/T, rho0/rho
%   Written: Chris Kniffin, Spring 2016
g = 1.4;
p0p = (1+(g-1)/2*M^2)^(g/(g-1));
A_As = 1/M*(2/(g+1)*(1+(g-1)/2*M^2))^((g+1)/(2*(g-1)));
T0T = 1+(g-1)/2*M^2;
rho0rho = (1+(g-1)/2*M^2)^(1/(g-1));

end

