% SOLUTION FOR THE PROJECT IN FEM 2020
% by Mika Persson & Wilhelm Treschow

%Problem b) (run 'problem_a.m' before running this file)
%% GENERATE C-MATRIX

C = zeros(ndof);
for elnr=1:nelm
    subdomain = t(4,elnr);         % what subdomain the current element belongs to
    x;
    switch subdomain
    case 1
        disp('negative one')
    case 2
        disp('zero')
    case 3
        disp('positive one')
    case 4
        
    end
end


