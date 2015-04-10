%% Creates the GCL parametric surfaces with 2 million points

function [x,y,z] = layer_eq_GCL_2(layer)
% Define rotation angles
xdeg= -35;
ydeg= 0;
zdeg= 0;

% Create mesh grid
ulin=linspace(pi*1/100,pi*98/100,2000);
vlin=linspace(pi*-23/100,pi*142.5/100,1000);
[u, v]=meshgrid(ulin,vlin);

% Define equations
x=(5.3-1*sin(u)+(1.00+layer*0.138)*cos(v)).*cos(u)*-500;
y=((5.5-2*sin(u)+(0.9+layer*0.114)*cos(v)).*sin(u))*750;
z=(sin(v+-0.13*(pi-u)))*(663+layer*114)+2500*sin(u);

% Create arrays and rotate
M=[reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)];
[x,y,z] = rotate3(M,xdeg,ydeg,zdeg);
