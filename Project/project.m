% SOLUTION FOR THE PROJECT IN FEM 2020

%% PREPROCESSOR
clear
clc

mesh = load('mesh.mat');  % geometry/mesh of half the IC
p = mesh.p;  % points/nodes (x:y for each column)
e = mesh.e;  % edges  (rows 1,2: node # of el. seg., 5: edge label bd. seg.)
t = mesh.t;  % triangles (4xnelm, rows 1-3: node numbers, 4: subdomain)

% DONT FORGET TO CONVERT FROM mm TO m!!

domain_k = [385; 149; 5; 385];  % k for each domain domain
thickness = 10;                 % thickness of the IC
initTemp = 30;                  % initial temperature
envTemp = 18;                   % environment temperature
Q = 5*10^7;                     % heat generated by the die
domain_Q = [0;Q;0;0];           % Q in the four domains

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

[ex,ey]=coordxtr(edof,coord,dof,3);  % x- and y coordinates for each element

% Check which segments that should have convections
er = e([1 2 5],:);        % reduced e (only interested of rows 1, 2 and 5)
conv_segments = [1 9 18]; % choosen boundary segments
edges_conv = [];          % edges with convection

for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];  % may use preallocation
    end
end

%% SOLVER
ndof = 2*nnod;       % number degrees of freedom (x & y per node)
D = eye(2);          % 2D constitutive element matrix
K = zeros(ndof);     % stiffness matrix
f = zeros(ndof, 1);  % force vector 

for elnr = 1:nelm
    
    subdomain = t(4,elnr);
    temp_k = domain_k(subdomain);
    temp_Q = domain_Q(subdomain);
    temp_D = temp_k * D;
    
    [Ke, fe] = flw2te(ex(elnr,:), ey(elnr,:), thickness, temp_D, temp_Q);
    indx = edof(elnr,2:end);              
    disp(indx)
    K(indx,indx) = K(indx,indx)+Ke; 
    
    indx = edof(elnr,2:end);               
    f(indx) = f(indx) + fe;             
end

convection_nodes = unique(edges_conv);  % nodes that belong to the convection boundary
nr_rows = size(convection_nodes, 1);    % number of nodes on the convection boundary
init_temp = zeros(nr_rows, 1);          % initial temperature on conv. bound.

bc = [convection_nodes, init_temp];

T = solveq(K, f, bc);  % nodal temperatures

%% POST PROCESSOR
eT=extract(edof,T);   % element temperatures

% In order to plot both sides of the symmetry cut:
patch(ex',ey',eT');   %
hold on     
patch(-ex',ey',eT');  % 

title('Temperature distribution [C]')
colormap(jet);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal

%patch(ex�,ey�,eT�,�EdgeColor�,�none�)  % removes the mesh lines

% SOME MORE PLOTTING STUFF THAT I HAVEN'T COPIED FROM THE PDF



