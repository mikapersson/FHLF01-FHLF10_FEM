close all
clear all

% SOLVE 1-D HEAT FLOW PROBLEM
A=10; Q = 100; k=5; q = 15;   % given in the problem
xStart = 2; xEnd = 8;         % where the object starts and ends
nelm=500; nen=2; ndof= nelm+1;  % #elements, #endpoints/element, #degoffreed.

x = linspace(xStart, xEnd, ndof);  % divide the object
L = (x(end)-x(1))/nelm;            % length of each element

Edof = zeros(nelm, nen+1);

for elnr=1:nelm
   Edof(elnr,:) = [elnr elnr elnr+1]; 
end

Ke = spring1e(A*k/L);                % element stiffness matrix
K = zeros(ndof);                     % stiffness matrix
Fb = zeros(ndof,1); Fb(end) = -q*A;  % boundary vector
Fl = zeros(ndof, 1);                 % load vector

for i=1:nelm
    Fl(i) = Fl(i) + Q*L/2;
    Fl(i+1) = Fl(i+1) + Q*L/2;
end

F = Fb + Fl;

for elnr = 1:nelm
    K = assem(Edof(elnr, :), K, Ke);
end

bc = [1 0];
T = solveq(K,F,bc);

plot(x,T)

% ANSWERS:
% 3 elements: T = [0 14 20 18]'
% 20 elements: T = [0;2.61000000000000;5.04000000000000;7.29000000000001;
% 9.36000000000001;11.2500000000000;12.9600000000000;14.4900000000000;
% 15.8400000000000;17.0100000000000;18.0000000000000;18.8100000000000;
% 19.4400000000000;19.8900000000000;20.1600000000001;20.2500000000001;
% 20.1600000000000;19.8900000000000;19.4400000000000;18.8100000000000;
% 18.0000000000000]
