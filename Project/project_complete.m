% SOLUTION FOR THE PROJECT IN FEM 2020
% by Mika Persson & Wilhelm Treschow

% Contains solutions to all subproblems (see sections)

% Problem a)
%% PREPROCESSOR 
clear
clc

mesh = load('mesh_fine2.mat');  % geometry/mesh of half the IC
p = mesh.p;  % points/nodes (x:y for each column)
e = mesh.e;  % edges (rows 1,2: node # of el. seg., 5: edge label bd. seg.)
t = mesh.t;  % triangles (4xnelm, rows 1-3: node numbers, 4: subdomain)

% Constants
domain_k = [385; 149; 5; 385];  % k for each domain domain (domains are given in the geometry)
thickness = 10*10^-3;           % thickness of the IC
init_Temp = 30;                 % initial temperature
env_Temp = 18;                  % environment temperature
alpha_c = 40;                   % convection parameter 
Q = 0.75*5*10^7;                     % heat generated by the die (mult. w. 0.75 for effect of COVID-19)
domain_Q = [0;Q;0;0];           % Q in the four domains (domains are given in the geometry)

% Material constants
E_Ag = 7;
E_Si = 165;
E_Cu = 128;
nu_Ag = 0.3;
nu_Si = 0.22;
nu_Cu = 0.36;
alpha_Ag = 4*10^-5;
alpha_Si = 2.6*10^-6;
alpha_Cu = 17.6*10^-6;
rho_Ag = 2500;
rho_Si = 2530;
rho_Cu = 8930;
c_Ag = 1000;
c_Si = 703;
c_Cu = 386;
k_Ag = 5;
k_Si = 149;
k_Cu = 385;

% Calculate edof matrices
coord = p';
coord = coord*10^-3;                % scale from mm to m
enod=t(1:3,:)';                     % nodes of elements
nelm=size(enod,1);                  % number of elements
nnod=size(coord,1);                 % number of nodes
dof=(1:nnod)';                      % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];  % edof matrix for elasticity problem
    edof(ie,:)=[ie,enod(ie,:)];  
end

[ex,ey]=coordxtr(edof,coord,dof,3);  % x- and y coordinates for each node

er = e([1 2 5],:);        % reduced e (only interested of rows 1, 2 and 5)

% Check which segments that should have convections
conv_segments = [1 9 18]; % choosen boundary segments
edges_conv = [];          % edges with convection

% Construct/Fill 'edges_conv'
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];  
    end
end

nr_conv_edges = size(edges_conv,2);
nodes_conv = unique(edges_conv);             % nodes that belong to the convection boundary
nr_conv_nodes = size(nodes_conv, 1);         % number of nodes on the convection boundary
edges_conv_coord = zeros(nr_conv_nodes, 4);  % (x,y)-coordinate for the two nodes
                                             % for every edge represented
                                             % as [x1, x2, y1, y2]
% Fill edges_conv_coord-matrix
for edge_nr = 1:nr_conv_edges
    edges_conv_coord(1,edge_nr) = coord(edges_conv(1,edge_nr),1);  %x1
    edges_conv_coord(2,edge_nr) = coord(edges_conv(2,edge_nr),1);  %x2
    edges_conv_coord(3,edge_nr) = coord(edges_conv(1,edge_nr),2);  %y1
    edges_conv_coord(4,edge_nr) = coord(edges_conv(2,edge_nr),2);  %y2
end

%% SOLVER
ndof = nnod;         % number degrees of freedom (x & y per node)
D = eye(2);          % 2D constitutive element matrix
K = zeros(ndof);     % stiffness matrix
f = zeros(ndof, 1);  % force vector 

% Create the stiffness matrix (K) and force vector (f)
for elnr = 1:nelm
    
    subdomain = t(4,elnr);         % what subdomain the current element belongs to
    temp_k = domain_k(subdomain);  % current k-value
    temp_Q = domain_Q(subdomain);  % current Q-value
    temp_D = temp_k*D;             % current D-matrix
    
    % Get element stiffness matrix and force vector
    [Ke, fe] = flw2te(ex(elnr,:), ey(elnr,:), thickness, temp_D, temp_Q);
    
    if subdomain == 3  % only subdomain 3 has convection
        % Examine if the current element (elnr) has any boundaries with
        % convection, if so we shall modify Ke
        for i = 1:size(edges_conv,2)
            if ismember([edges_conv(1,i), edges_conv(2,i)],enod(elnr,:))
                 x1 = edges_conv_coord(1,i);
                 x2 = edges_conv_coord(2,i);
                 y1 = edges_conv_coord(3,i);
                 y2 = edges_conv_coord(4,i);
                 boundaryLength = sqrt((x2-x1)^2+(y2-y1)^2);
                 a = thickness*alpha_c*[0, 0, 0; 0, boundaryLength/3, boundaryLength/6; 0, boundaryLength/6, boundaryLength/3];
                 Ke = Ke+a;
                 fe = fe+thickness*alpha_c*env_Temp*[0; boundaryLength/2; boundaryLength/2];
            end
        end
        
    end
    
    % Assembling matrices
    index = edof(elnr,2:end);       
    K(index,index) = K(index,index)+Ke;              
    f(index) = f(index) + fe;     
end

T = solveq(K, f);  % nodal temperatures

%% POST PROCESSOR
eT=extract(edof,T);   % element temperatures

% Print maximal temperature:
% 68.2 C without COVID
% 55.7 C with COVID
disp(['Maximal Temperature: ', num2str(max(max(eT)))])

% In order to plot both sides of the symmetry cut:
patch(ex',ey',eT','EdgeColor','none');   
hold on     
patch(-ex',ey',eT','EdgeColor','none');  

title('Temperature distribution with reduced video quality [C]')
grid on
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal


%% Problem b) (run the sections above before running this section)
% Generate C-matrix
C = zeros(ndof);
for elnr=1:nelm
    subdomain = t(4,elnr);  % what subdomain the current element belongs to
    x = null(1);
    switch subdomain
    case 1
        x = c_Cu*rho_Cu;
    case 2
        x = c_Si*rho_Si;
    case 3
        x = c_Ag*rho_Ag;
    case 4
        x = c_Cu*rho_Cu;
    end
    
    Ce = plantml(ex(elnr,:), ey(elnr,:), x);
    index = edof(elnr,2:end);   
    C(index,index) = C(index,index)+Ce;   
end

time_step = 1;  % step in time integration
end_time = 40*60;

transient_T = zeros(size(T,1),end_time);  % preallocate node temperatures for each time point
temp_T = init_Temp*ones(size(T));         % current node temperatures
transient_T(:,1) = temp_T;                % plug in initial values

% Calculate node temperatures for each time step
for time=2:end_time
    temp_T = (C+time_step*K)\(C*temp_T + time_step * f);
    transient_T(:,time) = temp_T;
end

%% Animate transient heat transfer
for time=1:end_time
    temp_T = transient_T(:,time);
    eT=extract(edof,temp_T);   % element temperatures

    % In order to plot both sides of the symmetry cut:
    patch(ex',ey',eT','EdgeColor','none');   %
    hold on     
    patch(-ex',ey',eT','EdgeColor','none');  % 

    title('Temperature distribution [C]')
    colormap(hot);
    colorbar;
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    axis equal
    
    pause(0.01)
end

%% Decide which time point you want to examine
time_point = 40*60;
eT=extract(edof,transient_T(:,time_point));   % element temperatures

% Plot temperature at time 'time_point'
patch(ex',ey',eT','EdgeColor','none');   
hold on     
patch(-ex',ey',eT','EdgeColor','none');   

title('Temperature distribution at initial time point [C]')
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
    

%% Problem c) (run the sections above before running this section)
%
ndof = 2*nnod;        % number of degrees of freedom
E_values = 10^9 * [E_Cu E_Si E_Ag E_Cu];
v_values = [nu_Cu nu_Si nu_Ag nu_Cu];
alfa_values = [alpha_Cu alpha_Si alpha_Ag alpha_Cu];

ptype = 2;               % because of plain strain problem
Ks = zeros(ndof);        % elasticity stiffness matrix
fs = zeros(ndof, 1);     % elasticity force vector
Ep = [ptype thickness];  % from CALFEM manual p.87

element_DeltaT = zeros(nelm,1);  % mean temperature difference for every element

% Create Ks and fs (and element_DeltaT)
for i =1:nelm
   
    % The temperature difference in each node
    DeltaT1 = T(t(1,i))-init_Temp;             
    DeltaT2 = T(t(2,i))-init_Temp;
    DeltaT3 = T(t(3,i))-init_Temp;
    DeltaT = (DeltaT1 + DeltaT2 + DeltaT3)/3;  % this is the mean temperatuer difference for an element

    element_DeltaT(i,:) = DeltaT;
    
    subdomain = t(4,i);  
    
    % Get constants for current element
    E = E_values(subdomain);  
    v = v_values(subdomain);
    alfa = alfa_values(subdomain);
    
    Ds = E/((1+v)*(1-2*v))*[1-v v 0; v 1-v 0; 0 0 0.5*(1-2*v)];  % see p.255
    Ds0 = Ds*(1+v)*alfa*DeltaT*[1; 1; 0];  % D*epsilon_0 (thermal strain)
    
    Kse=plante(ex(i,:),ey(i,:),Ep,Ds);     % element Ks
    fse=plantf(ex(i,:),ey(i,:),Ep,Ds0');   % element fs
    
    % Assembling matrices
    indx = edof_S(i,2:end);
    Ks(indx,indx) = Ks(indx,indx)+Kse;
    fs(indx) = fs(indx) + fse;
    
end

symmetry_edge_labels = [10,11,12,13];  % boundary segments at symmetry axis (x=0)
fixed_edge_label = 7;  % fixated bottom boundary segment (y=0)
bc = [];               
found_nodes = [];      % keeps track of visited nodes

% Construct boundary value vector (bc) (not so pretty..)
for ei=1:size(er,2)
   node_1 = er(1,ei); 
   node_2 = er(2,ei); 
   edge_label = er(3,ei); 
    if sum(ismember(symmetry_edge_labels, edge_label)) ~= 0  % if the edge is on the symmetry axis
        if sum(ismember(found_nodes, node_1)) == 0           % we haven't added node1 to bc
           new_bc = [dof_S(node_1,1), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_1];
        end
        if sum(ismember(found_nodes, node_2)) == 0           % we haven't added node2 to bc
           new_bc = [dof_S(node_2,1), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_2];
        end
        
    elseif edge_label == fixed_edge_label                    % if we are at the bottom (y=0)
        if sum(ismember(found_nodes, node_1)) == 0           % we haven't added node1 to bc
           new_bc = [dof_S(node_1,1), 0; dof_S(node_1,2), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_1];
        end
        if sum(ismember(found_nodes, node_2)) == 0          % we haven't added node2 to bc
           new_bc = [dof_S(node_2,1), 0; dof_S(node_2,2), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_2];
        end
    end
end

u = solveq(Ks,fs,bc);    % solve for 2D-displacement (at each node)
ed = extract(edof_S,u);  % element displacements

%% Plot displacement
% Calculate displaced coordinates 
mag = 100;  % Magnification (due to small deformations)
exd = ex + mag*ed(:,1:2:end);
eyd = ey + mag*ed(:,2:2:end);

% Plot deformed IC over original circuit
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3) 
hold on 
patch(-ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3) 
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3)
hold on
patch(-exd',eyd',[0 0 0],'FaceAlpha',0.3)

xlabel('x-position [m]')
ylabel('y-position [m]')

axis equal

title('Displacement field with reduced video quality [Magnitude enhancement 100]')

%% Von mises stress-field (similar to the computations as above)

Seff_el = zeros(nelm,1);
for j = 1:nelm
    subdomain = t(4,j);
    
    E = E_values(subdomain);
    v = v_values(subdomain);
    alfa = alfa_values(subdomain);
    DeltaT = element_DeltaT(j,:);
    
    Ds = E/((1+v)*(1-2*v))*[1-v v 0; v 1-v 0; 0 0 0.5*(1-2*v)];
    
    [es,et] = plants(ex(j,:),ey(j,:),Ep, Ds, ed(j,:));
    
    Ds0 = Ds*(1+v)*alfa*DeltaT*[1; 1; 0;]; 
    es = es - transpose(Ds0);  % correction for thermal stress
    
    sigma_xx = es(1);  
    sigma_yy = es(2);  
    tau_xy = es(3);
    
    
    epsilon_xx = et(1);
    epsilon_yy = et(2);
    gamma_xy = et(3);
    
    sigma_zz = v*(sigma_xx + sigma_yy) - alfa*DeltaT*E;
    
    % This is the von Mises stress field for an element
    Seff_el(j) = sqrt(sigma_xx^2 + sigma_yy^2 + sigma_zz^2 - sigma_xx*sigma_yy - sigma_xx*sigma_zz - sigma_yy*sigma_zz + 3*tau_xy^2);
    
end

%% Plot stress field

Seff_nod = zeros(nnod,1);
for i=1:size(coord,1) 
    [c0,c1]=find(edof(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
end

Seff = extract(edof,Seff_nod);  % element stresses

% Get peak stress
peak_stress = max(Seff_el);
disp(['Peak Stress: ', num2str(peak_stress)])
% 146 MPa without reduction
% 97.9 MOa with reduction

% Get element with peak stress
el_nr = find(Seff_el == peak_stress);
disp(['at element: ', num2str(el_nr)])
% 43

% Get coordinates of peak stress
point_nr = t(1,el_nr);
x_coord = p(1,point_nr);
y_coord = p(2,point_nr);
disp(['at coordinates: ', num2str(x_coord), ':', num2str(y_coord)])
% x=3.5, y=0


% Plot stress field of deformed IC
patch(exd',eyd',Seff','EdgeColor','none');   
hold on     
patch(-exd',eyd',Seff','EdgeColor','none');  

title('Stress field with reduced video quality [Pa]')
colormap(jet);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal


















