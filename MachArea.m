function [ Msub,Msup ] = MachArea( AAs )
%Input: A/Astar Outputs: Sub- and super-sonic machnumbers
%   Written: C. Kniffin
Msubl = 0; % lower bound of sub sonic mach number
Msubu = .999999999; % upper bound of sub sonic mach number
Msupl = 1.00000001; % lower bound of super sonic mach number
Msupu = 20; % upper bound of super sonic mach number
Msub = Msubl;
Msup = Msupl;
g = 1.4;

AAseqn = 0;
k = 0;
kmax = 100;

while abs(AAseqn - AAs^2) > eps && k < kmax
    k = k+1;
    AAseqn = 1/Msub^2*(2/(g+1)*(1+(g-1)/2*Msub^2))^((g+1)/(g-1));
    if AAseqn < AAs^2
        Msubu = Msub;
        Msub  = (Msubl+Msub)/2;
    else
        Msubl = Msub;
        Msub  = (Msubu+Msub)/2;
    end
end

AAseqn = 0;
k = 0;
kmax = 100;

while abs(AAseqn - AAs^2) > eps && k < kmax
    k = k+1;
    AAseqn = 1/Msup^2*(2/(g+1)*(1+(g-1)/2*Msup^2))^((g+1)/(g-1));
    if AAseqn > AAs^2
        Msupu = Msup;
        Msup  = (Msupl+Msup)/2;
    else
        Msupl = Msup;
        Msup  = (Msupu+Msup)/2;
    end
end