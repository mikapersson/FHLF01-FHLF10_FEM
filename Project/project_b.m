% SOLUTION FOR THE PROJECT IN FEM 2020
% by Mika Persson & Wilhelm Treschow

%Problem b) (run the sections above before running this section)

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
end_time = 1000;

transient_T = zeros(size(T,1),end_time);  % node temperatures for each time point
temp_T = init_Temp*ones(size(T));
transient_T(:,1) = temp_T;  % node temperatures for each time point

% Calculate node temperatures for each time step
for t=2:end_time
    temp_T = (C+time_step*K)\(C*temp_T + time_step * f);
    transient_T(:,t) = temp_T;
end


% Animate transient heat transfer
for t=1:end_time
    temp_T = transient_T(:,t);
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
time_point = 1000;
eT=extract(edof,transient_T(:,time_point));   % element temperatures

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
    









