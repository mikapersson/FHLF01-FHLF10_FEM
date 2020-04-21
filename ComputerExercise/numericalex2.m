
% SOLVE 2-D HEAT FLOW PROBLEM

clf

% defined coordinates

nf=7;
nr=6;  % must be even
R1=2;
R2=10;
coord=[];
for R=R1:(R2-R1)/nf:R2
for fi=0:pi/2/nr:pi/2
   if fi<pi/4
     Rf=(R1-R)/(R1-R2)*R/cos(fi)+(R2-R)/(R2-R1)*R;
     coord=[coord;Rf*cos(fi) Rf*sin(fi)];
   else
     Rf=(R1-R)/(R1-R2)*R/sin(fi)+(R2-R)/(R2-R1)*R;
     coord=[coord;Rf*cos(fi) Rf*sin(fi)];
   end
end
   
end
plot(coord(:,1),coord(:,2),'*')
axis equal
hold on

% global dof for nodal points

ndof=0;
dof=[];
for R=R1:(R2-R1)/nf:R2
  for fi=0:pi/2/nr:pi/2
    ndof=ndof+1;
    dof=[dof; ndof];
  end
end

% topology matriz

nelm=0;

edof=[];
for fn=0:nf-1
 
  for ir=0:nr-1
     n1=fn*(nr+1)+1+ir;
     n2=fn*(nr+1)+2+ir;
     n3=(fn+1)*(nr+1)+1+ir;
     n4=(fn+1)*(nr+1)+2+ir;    

     nelm=nelm+1;
     edof=[edof; nelm n3 n2 n1];
     nelm=nelm+1;
     edof=[edof; nelm n3 n4 n2];
  end
end

[ex,ey]=coordxtr(edof,coord,dof,3);
eldraw2(ex,ey,[1 2 0],edof(:,1))
%pause

% Define boundary conditions
%
% find nodal points on boundaries

bc3=find(abs(coord(:,1)-R2)<1e-3);
bc4=find(abs(coord(:,2)-R2)<1e-3);
bcarc=find(abs(coord(:,1).^2+coord(:,2).^2-R1^2)<1e-3);


bcarc(:,2)=1000; % circle
bc3(:,2)=100;    % top boundary
bc4(:,2)=100;    % right boundary

bc=[bcarc;bc3;bc4];


bc(11,:)=[]; % Remove 