% This script is written and read by pdetool and should NOT be edited.
% There are two recommended alternatives:
 % 1) Export the required variables from pdetool and create a MATLAB script
 %    to perform operations on these.
 % 2) Define the problem completely using a MATLAB script. See
 %    http://www.mathworks.com/help/pde/examples/index.html for examples
 %    of this approach.
function pdemodel
[pde_fig,ax]=pdeinit;
pdetool('appl_cb',1);
set(ax,'DataAspectRatio',[4 5.25 1]);
set(ax,'PlotBoxAspectRatio',[0.92736077481840218 1 0.92736077481840218]);
set(ax,'XLim',[-2 6]);
set(ax,'YLim',[-0.5 10]);
set(ax,'XTick',[ -2,...
 -1.5,...
 -1,...
 -0.5,...
 0,...
 0.5,...
 1,...
 1.5,...
 2,...
 2.5,...
 3,...
 3.5,...
 4,...
 4.5,...
 5,...
 5.5,...
 6,...
 6.5,...
 7,...
 7.5,...
 8,...
]);
set(ax,'YTick',[ 0,...
 0.5,...
 1,...
 1.5,...
 2,...
 2.5,...
 3,...
 3.5,...
 4,...
 4.5,...
 5,...
 5.5,...
 6,...
 6.5,...
 7,...
 7.5,...
 8,...
 8.5,...
 9,...
 9.5,...
 10,...
]);
pdetool('gridon','on');

% Geometry description:
pderect([0 4 7 3],'R1');
pdepoly([ 3,...
 4,...
 4,...
],...
[ 7,...
 7,...
 6,...
],...
 'P1');
pderect([0 1 4.5 4],'R2');
pderect([0 0.5 5 4.5],'R3');
pdepoly([ 2,...
 2,...
 3.25,...
 3.5,...
 3.5,...
 4,...
 4,...
 3.5,...
],...
[ 4.5,...
 4,...
 4,...
 3.75,...
 0,...
 0,...
 4,...
 4.5,...
],...
 'P2');
set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','(R1-P1)+R2+R3+P2')

% Boundary conditions:
pdetool('changemode',0)
pdetool('removeb',[22 ]);
pdesetbd(24,...
'dir',...
1,...
'1',...
'0')
pdesetbd(22,...
'dir',...
1,...
'1',...
'0')
pdesetbd(18,...
'dir',...
1,...
'1',...
'0')
pdesetbd(17,...
'dir',...
1,...
'1',...
'0')
pdesetbd(16,...
'dir',...
1,...
'1',...
'0')
pdesetbd(13,...
'dir',...
1,...
'1',...
'0')
pdesetbd(12,...
'dir',...
1,...
'1',...
'0')
pdesetbd(11,...
'dir',...
1,...
'1',...
'0')
pdesetbd(10,...
'dir',...
1,...
'1',...
'0')
pdesetbd(9,...
'dir',...
1,...
'1',...
'0')
pdesetbd(7,...
'dir',...
1,...
'1',...
'0')
pdesetbd(1,...
'dir',...
1,...
'1',...
'0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
setappdata(pde_fig,'MesherVersion','preR2013a');
pdetool('initmesh')

% PDE coefficients:
pdeseteq(1,...
'1.0',...
'0.0',...
'10.0',...
'1.0',...
'0:10',...
'0.0',...
'0.0',...
'[0 100]')
setappdata(pde_fig,'currparam',...
['1.0 ';...
'0.0 ';...
'10.0';...
'1.0 '])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
char('0','1000','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 0 0 0 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')
