% SOLUTION FOR THE PROJECT IN FEM 2020

clear
clc

mesh = load('mesh.mat');  % geometry/mesh of half the IC
p = mesh.p;  % points/nodes (x:y for each column)
e = mesh.e;  % edges  (rows 1,2: node # of el. seg., 5: edge label bd. seg.)
t = mesh.t;  % triangles (4xnelm, rows 1-3: node numbers, 4: subdomain)

coord = p';
enod=t(1:3,:)';                     % nodes of elements
nelm=size(enod,1);                  % number of elements
nnod=size(coord,1);                 % number of nodes
dof=(1:nnod)';                      % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number
for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end