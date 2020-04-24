% SOLUTION FOR THE PROJECT IN FEM 2020

%% PREPROCESSOR
clear
clc

mesh = load('mesh.mat');  % geometry/mesh of half the IC
p = mesh.p;  % points/nodes (x:y for each column)
e = mesh.e;  % edges  (rows 1,2: node # of el. seg., 5: edge label bd. seg.)
t = mesh.t;  % triangles (4xnelm, rows 1-3: node numbers, 4: subdomain)

thickness = 10;  % thickness of the IC
initTemp = 30;   % initial temperature
envTemp = 18;    % environment temperature
Q = 5*10^7;      % heat generated by the die

% Calculate edof matrices
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

% Check which segments that should have convections
er = e([1 2 5],:); % Reduced e (only interested of rows 1, 2 and 5)
conv_segments = [1 9 17]; % Choosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];  % may use preallocation
    end
end

%% SOLVER
Kte = 0;                            % element stiffness matrix (triangle element)
Kt = assem(edof(el,:),Kt,Kte);      % stiffness matrix
indx = edof(el,2:end);              % 
Kt(indx,indx) = Kt(indx,indx)+Kte;  % 
f = insert(edof(el,:),f,fe);        % 
indx = edof(el,2:end);              % 
f(indx) = f(indx) + fe;             % 

%T = soleq...                 % nodal temperatures

%% POST PROCESSOR
eT=extract(edof,T);   % element temperatures

% In order to plot both sides of the symmetry cut:
patch(ex',ey',eT');   %
hold on     
patch(-ex',ey',eT');  % 

%patch(ex�,ey�,eT�,�EdgeColor�,�none�)  % removes the mesh lines

% SOME MORE PLOTTING STUFF THAT I HAVEN'T COPIED FROM THE PDF


