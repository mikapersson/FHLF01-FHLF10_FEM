% SOLUTION FOR THE PROJECT IN FEM 2020
% by Mika Persson & Wilhelm Treschow

% Problem c) (run the sections above before running this section)
%%
ndof = 2*nnod;                        % number of degrees of freedom
E_values = 10^9 * [E_Cu E_Si E_Ag E_Cu];
v_values = [nu_Cu nu_Si nu_Ag nu_Cu];
alfa_values = [alpha_Cu alpha_Si alpha_Ag alpha_Cu];

ptype = 2;            % because of plain strain problem
Ks = zeros(ndof);     % elasticity stiffness matrix
fs = zeros(ndof, 1);  % elasticity force vector
Ep = [ptype ep];      % from CALFEM manual p.87

element_DeltaT = zeros(nelm,1);  

% Create Ks and fs (and element_DeltaT)
for i =1:nelm
   
    DeltaT1 = T(t(1,i))-init_Temp;             % the temperature difference in each node
    DeltaT2 = T(t(2,i))-init_Temp;
    DeltaT3 = T(t(3,i))-init_Temp;
    DeltaT = (DeltaT1 + DeltaT2 + DeltaT3)/3;  % this is the mean temperatuer difference for an element

    element_DeltaT(i,:) = DeltaT;
    
    subdomain = t(4,i);  
    
    % get constants for current element
    E = E_values(subdomain);  
    v = v_values(subdomain);
    alfa = alfa_values(subdomain);
    
    Ds = E/((1+v)*(1-2*v))*[1-v v 0; v 1-v 0; 0 0 0.5*(1-2*v)];  % see p.255
    Ds0 = Ds*(1+v)*alfa*DeltaT*[1; 1; 0];
    
    Kse=plante(ex(i,:),ey(i,:),Ep,Ds);    % element Ks
    fse=plantf(ex(i,:),ey(i,:),Ep,Ds0');  % element fs
    
    indx = edof_S(i,2:end);
    Ks(indx,indx) = Ks(indx,indx)+Kse;
    fs(indx) = fs(indx) + fse;
    
end

symmetry_edge_labels = [10,11,12,13];
fixed_edge_label = 7;
bc = [];
found_nodes = [];  % keeps track of visited nodes

% Construct boundary value vector (bc)  (not so pretty..)
for ei=1:size(er,2)
   [node_1, node_2, edge_label] = er(:,ei); 
    if ismember(symmetry_edge_labels, edge_label)  % if the edge is on the symmetry axis
        if ismember(found_nodes, node_1)           % we haven't added node1 to bc
           new_bc = [dof_S(node_1,1), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_1];
        end
        if ismember(found_nodes, node_2)           % we haven't added node2 to bc
           new_bc = [dof_S(node_2,1), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_2];
        end
        
    elseif edge_label == fixed_edge_label          % if we are at the bottom
        if ismember(found_nodes, node_1)           % we haven't added node1 to bc
           new_bc = [dof_S(node_1,1), 0; dof_s(node_1,2), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_1];
        end
        if ismember(found_nodes, node_2)           % we haven't added node2 to bc
           new_bc = [dof_S(node_2,1), 0; dof_S(node_2,2), 0];
           bc = [bc; new_bc];
           found_nodes = [found_nodes, node_2];
        end
    end
end

u = solveq(Ks,fs,bc);
ed = extract(edof_S,u);


%% Plot displacement
% Calculate displaced coordinates 
mag = 100;              % Magnification (due to small deformations)
exd = ex + mag*ed(:,1:2:end);
eyd = ey + mag*ed(:,2:2:end);

figure() 
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3) 
hold on 
patch(-ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3) 
hold on
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3)
hold on
patch(-exd',eyd',[0 0 0],'FaceAlpha',0.3)

axis equal

title('Displacement field [Magnitude enhancement 100]')

%% Von mises stress-field

stress = zeros(nelm,4);
vonmises_field = zeros(nelm,1);
for j = 1:nelm
    subdomain = t(4,i);
    
    E = E_values(subdomain);
    v = v_values(subdomain);
    alfa = alfa_values(subdomain);
    DeltaT = element_DeltaT(j,:);
    
    Ds = E/((1+v)*(1-2*v))*[1-v v 0; v 1-v 0; 0 0 0.5*(1-2*v)];
    
    [es,et] = plants(ex(j,:),ey(j,:),Ep, Ds, ed(j,:));
    
    %plane stress: epsilon_zz = 0
    
    %es = [sigma_xx sigma_yy sigma_xy]
    
    %et = [epsilon_xx epsilon_yy gamma_xy]
    
    Ds0 = Ds*(1+v)*alfa*DeltaT*[1; 1; 0;]; 
    es = es - transpose(Ds0);
    
    sigma_xx = es(1);
    sigma_yy = es(2);
    tau_xy = es(3);
    
    
    epsilon_xx = et(1);
    epsilon_yy = et(2);
    gamma_xy = et(3);
    
    %sigma_zz = E*v/((1+v)*(1-2*v))*(epsilon_xx + epsilon_yy) - alfa*E*DeltaT/(1-2*v);
    sigma_zz = v*(sigma_xx + sigma_yy) - alfa*DeltaT*E;
    
    %sigma_zz = E*v/((1+v)*(1-2*v))*(et(1) + et(2)) - alfa*E*DeltaT/(1-2*v)
    
    
    %stress(j,:) = stress(j,:) + [es-Ds0 sigma_zz];
    
    % This is the von Mises stress field for an element
    %vonmises_field(j) = sqrt(stress(j,1)^2 + stress(j,2)^2 + stress(j,3)^2 - stress(j,1)*stress(j,2) - stress(j,1)*stress(j,3) - stress(j,2)*stress(j,3) + 3*stress(j,4)^2);
    vonmises_field(j) = sqrt(sigma_xx^2 + sigma_yy^2 + sigma_zz^2 - sigma_xx*sigma_yy - sigma_xx*sigma_zz - sigma_yy*sigma_zz + 3*tau_xy^2);
    
end

%% Plot stress field

Seff_el = vonmises_field;

Seff_nod = zeros(nnod,1);
for i=1:size(coord,1) 
    [c0,c1]=find(edof(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
end

Seff = extract(edof,Seff_nod);
patch(ex',ey',Seff','EdgeColor','none');   %
hold on     
patch(-ex',ey',Seff','EdgeColor','none');  % 

title('Stress field [Pa]')
colormap(jet);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal












