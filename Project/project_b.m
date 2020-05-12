% SOLUTION FOR THE PROJECT IN FEM 2020
% by Mika Persson & Wilhelm Treschow

%Problem b) (run 'problem_a.m' before running this file)
%% GENERATE C-MATRIX

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

time_step = 1;  % step i time integration

transient_T = zeros(size(T,1),20);  % node temperatures for each time point
temp_T = init_Temp*ones(size(T));
transient_T(:,1) = temp_T;  % node temperatures for each time point

% Calculate node temperatures for each time step

temp_f = 
for t=2:20
    temp_T = (C+time_step*K)\(C*temp_T + time_step * temp_f);
    trainsient_T(t,1) = temp_T;
end


eT=extract(edof,T);   % element temperatures

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
    









