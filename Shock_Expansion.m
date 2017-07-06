%% Calculates cl and cd using shock expansion theory %%
%  Requires: IsentropicInv.m, ExpAngle.m, ExpMach.m, Beta.m, Shock.m
%   Written: Chris Kniffin, Spring 2016
clc; clear; close all;

%% -------------- Geometry Creation ------------- %%
%  Coordinates and angles can be entered manually instead

% Number of Partitions
n = 10;

% Generate X-Coordinates
x = linspace(0,1,n+1);

% Generate Y-Coordinates
yl = 2.475 - sqrt(6.126-x.^2+x);
yu = -2.475 + sqrt(6.126-x.^2+x);

% Angle of Attack
alpha = 1.5;

% Calculate Angles wrt Freestream
for i = 1:n
    angu(i) = asind((yu(i+1)-yu(i))/(x(i+1)-x(i)))-alpha;
    angl(i) = asind((yl(i+1)-yl(i))/(x(i+1)-x(i)))-alpha;
end
diffangu = angu(1);
diffangl = angl(1);

% Calculate Angles wrt Local Flow
for i = 1:n-1
    diffangu(i+1) = angu(i+1)-angu(i);
    diffangl(i+1) = angl(i+1)-angl(i);
end

% Calculate Surface Lengths
for i = 1:n
    slu(i) = sqrt(abs(yu(i)-yu(i+1))^2+abs(x(i)-x(i+1))^2);
    sll(i) = sqrt(abs(yl(i)-yl(i+1))^2+abs(x(i)-x(i+1))^2);
end

% Total Length
totlenu = sum(slu);
totlenl = sum(sll);

%% ---------------------------------------------- %%
% Free Stream Conditions
Minf = 1.5;
g = 1.4;
Mu(1) = Minf;
Ml(1) = Minf;
p0u(1) = 1;
pu(1) = IsentropicInv( Mu );
p0l(1) = 1;
pl(1) = IsentropicInv( Ml );

% Calculate Upper and Lower Pressures
for i = 1:n
    if diffangu(i) < 0
        expru(i) = 1;
        p0u(i+1) = p0u(i);
        v1 = ExpAngle(Mu(i));
        v2 = v1 - diffangu(i);
        Mu(i+1) = ExpMach(v2);
        pp0 = IsentropicInv( Mu(i+1) );
        pu(i+1) = p0u(i)*pp0;
    else
        expru(i) = 0;
        B = Beta(Mu(i),diffangu(i));
        M1n = Mu(i)*sind(B);
        M1ns(i) = M1n;
        [ M2n,p2p1,p02p01] = Shock( M1n );
        Mu(i+1) = M2n/(sind(B-diffangu(i)));
        p0u(i+1) = p02p01*p0u(i);
        pu(i+1) = p2p1*pu(i);
    end
end
for i = 1:n
    if diffangl(i) > 0
        exprl(i) = 1;
        p0l(i+1) = p0l(i);
        v1 = ExpAngle(Ml(i));
        v2 = v1 + diffangl(i);
        Ml(i+1) = ExpMach(v2);
        pp0 = IsentropicInv( Ml(i+1) );
        pl(i+1) = p0l(i)*pp0;
    else
        exprl(i) = 0;
        B = Beta(Ml(i),-diffangl(i));
        M1n = Ml(i)*sind(B);
        [ M2n,p2p1,p02p01] = Shock( M1n );
        Ml(i+1) = M2n/(sind(B+diffangl(i)));
        p0l(i+1) = p02p01*p0l(i);
        pl(i+1) = p2p1*pl(i);
    end
end

% Total Drag and Lift Calculation
dl = 0;
du = 0;
ll = 0;
lu = 0;

for i = 1:n
    dl = dl - sll(i)*pl(i+1)*sind(angl(i));
    du = du + slu(i)*pu(i+1)*sind(angu(i));
    ll = ll + sll(i)*pl(i+1)*cosd(angl(i));
    lu = lu - slu(i)*pu(i+1)*cosd(angu(i));
end

% Normalize by length
l = ll/totlenl+lu/totlenu;
d = dl/totlenl+du/totlenu;

% Coefficients
cl = 2/(g*Minf^2)*l
cd = 2/(g*Minf^2)*d

% L/D ratio
ld = cl/cd