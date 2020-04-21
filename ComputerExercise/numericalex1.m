close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preprocessor
A=10; Q = 100;
k=1;
xStart = 2; xEnd = 8;
nelm=3; nen=2;
ndof=8;

Coord=[ -L 0  % coordinates for the nodes
-L/sqrt(2) L/sqrt(2)
0 0
L/sqrt(2) L/sqrt(2)];

Edof=[1 1 2 5 6
2 3 4 5 6
3 5 6 7 8];




Dof=[ 1 2
3 4
5 6
7 8];  

[Ex,Ey]=coordxtr(Edof,Coord,Dof,nen);

K=zeros(ndof);
F = zeros(ndof,1);

F(5)=10;
F(6)=1000000;

bc=[1 0
2 0
3 0
4 0
7 0
8 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Solver
for elnr=1:nelm
Ke=bar2e(Ex(elnr,:),Ey(elnr,:),[E A]);
K=assem(Edof(elnr,:),K,Ke);
end

a=solveq(K,F, bc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Post processor

ed=extract(Edof,a);
es=bar2s(Ex,Ey,[E A],ed);

eldraw2(Ex,Ey,[1 4 1])

eldisp2(Ex,Ey,ed,[1 2 1])