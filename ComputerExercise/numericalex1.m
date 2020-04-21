close all
clear all

% SOLVE 1-D HEAT FLOW PROBLEM
A=10; Q = 100; k=5; q = 15;   % given in the problem
xStart = 2; xEnd = 8;         % where the object starts and ends
nelm=2000; nen=2; ndof= nelm+1;  % #elements, #endpoints/element, #degoffreed.

x = linspace(xStart, xEnd, ndof);  % divide the object
L = (x(end)-x(1))/nelm;            % length of each element

Ke = A*k/L*[1 -1; -1 1];           % element stiffness matrix

Edof = zeros(nelm, nen+1);

for elnr=1:nelm
   Edof(elnr,:) = [elnr elnr elnr+1]; 
end

K=zeros(ndof);
F = zeros(ndof,1); F(end) = -q*A;
Fe = [Q*L/2 ; Q*L/2];

for elnr = 1:nelm
    [K, F] = assem(Edof(elnr, :), K, Ke, F, Fe);
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
