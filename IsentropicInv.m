function [ pp0,As_A,TT0,rhorho0 ] = IsentropicInv( M )
%Outputs: p/p0, As/A, T/T0, rho/rho0
%   Written: Chris Kniffin, Spring 2016
[ p0p,A_As,T0T,rho0rho ] = Isentropic( M );
pp0 = p0p^-1;
As_A = A_As^-1;
TT0 = T0T^-1;
rhorho0 = rho0rho^-1;
end

