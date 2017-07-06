function [ M ] = ExpMach( v )
%Outputs: M Inputs: v
%   Written: Chris Kniffin, Spring 2016

g = 1.4;
k = 0;
kmax = 1000;
Mmax = 90;
Mmin = 0;
M = Mmin;
veqn = 0;

while abs(v-veqn)>eps && k<kmax
    k = k+1;
    veqn = sqrt((g+1)/(g-1))*atand(sqrt(((g-1)/(g+1)*(M^2-1))))-...
        atand(sqrt(M^2-1));
    if veqn > v
        Mmax = M;
        M = (M + Mmin)/2;
    else
        Mmin = M;
        M = (M + Mmax)/2;
    end
end
end

