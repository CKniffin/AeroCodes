function [ M2,p2p1,p02p01,T2T1,rho2rho1 ] = Shock( M )
%Outputs: M2,p2/p1,p02/p01,T2/T1,rho2/rho1
%   Written: Chris Kniffin, Spring 2016
g = 1.4;
M2 = sqrt((1+(g-1)/2*M^2)/(g*M^2-(g-1)/2));
p2p1 = 1+2*1.4/2.4*(M^2-1);
p02p01 = (1+2*g/(g+1)*(M^2-1))^(-1/(g-1))*((g+1)*M^2/((g-1)*M^2+2))^(g/(g-1));
T2T1 = 1+2*(g-1)/(g+1)^2*(g*M^2+1)/M^2*(M^2-1);
rho2rho1 = (g+1)*M^2/((g-1)*M^2+2);
end

