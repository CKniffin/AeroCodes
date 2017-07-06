function [ B ] = Beta( M,theta )
%Outputs: Beta, Inputs: M, theta
%   Written: Chris Kniffin, Spring 2016

g = 1.4;
k = 0;
kmax = 1000;
Bmax = 67;
Bmin = 0;
B = Bmin;
tantheta = 0;

while abs(tand(theta)-tantheta)>eps && k<kmax
    k = k+1;
    tantheta = 2*cotd(B)*(M^2*sind(B)^2-1)/((M^2*(g+cosd(2*B))+2));
    if tantheta > tand(theta)
        Bmax = B;
        B = (B + Bmin)/2;
    else
        Bmin = B;
        B = (B + Bmax)/2;
    end

end

