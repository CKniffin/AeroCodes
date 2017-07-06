function [ v ] = ExpAngle( M )
%Outputs: v Inputs: M
%   Written: Chris Kniffin, Spring 2016
g = 1.4;

v = sqrt((g+1)/(g-1))*atand(sqrt(((g-1)/(g+1)*(M^2-1))))-...
    atand(sqrt(M^2-1));

end

