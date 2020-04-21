close all
clear all

% SOLVE 1-D HEAT FLOW PROBLEM
A=10; Q = 100; k=5; q = 15;   % given in the problem
xStart = 2; xEnd = 8;         % where the object starts and ends
nelm=3; nen=2; ndof= nelm+1;  % #elements, #endpoints/element, #degoffreed.

x = linspace(xStart, xEnd, ndof);  % divide the object
L = (x(end)-x(1))/nelm;            % length of each element

Ke = A*k/L*[1 -1; -1 1];           % element stiffness matrix

Edof = zeros(nelm, nen+1);

for elnr=1:nelm
   Edof(elnr,:) = [elnr elnr elnr+1]; 
end

K=zeros(ndof);
F = zeros(ndof,1); F(xEnd) = -q*A;
Fe = [Q*L/2 ; Q*L/2];

for elnr = 1:nelm
    [K, F] = assem(Edof(elnr, :), K, Ke, F, Fe);
end

bc = [1 0];
T = solveq(K,F,bc);
