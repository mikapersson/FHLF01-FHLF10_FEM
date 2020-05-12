% SOLUTION FOR THE PROJECT IN FEM 2020
% by Mika Persson & Wilhelm Treschow

%Problem b) (run 'problem_a.m' before running this file)
%% GENERATE C-MATRIX

C = zeros(ndof);
for elnr=1:nelm
    subdomain = t(4,elnr);         % what subdomain the current element belongs to
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



