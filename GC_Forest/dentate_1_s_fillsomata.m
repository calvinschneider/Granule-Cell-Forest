%% Create Grid of Packed Somata and Select Those That Lie Within GCL

function dentate_1_s_fillsomata()

% Define parameters
SomaSize    =   15.79;
Kepler      =   pi / (3 * sqrt (2));

% Create granule cell layer surface
[x_g1,y_g1,z_g1]    = layer_eq_GCL(-1.95);
[x_g2,y_g2,z_g2]    = layer_eq_GCL(-1.0);
[x_g3,y_g3,z_g3]    = layer_eq_GCL(0);
X_g                 = [x_g1;x_g2;x_g3];
Y_g                 = [y_g1;y_g2;y_g3];
Z_g                 = [z_g1;z_g2;z_g3];
[~,S_g]             = alphavol([X_g(:),Y_g(:),Z_g(:)],120);

% Define limits for granule cell somata grid
xmin    =   -3200;
xmax    =   3200;
ymin    =   -500;
ymax    =   4200;
zmin    =   -600;
zmax    =   1300;

% Round limits to nearest multiple of soma size
xmax2    = round((xmax-xmin)/SomaSize)*SomaSize + xmin;
ymax2    = round((ymax-ymin)/SomaSize)*SomaSize + ymin;
zmax2    = round((zmax-zmin)/(Kepler*SomaSize))*(Kepler*SomaSize) + zmin;
zlayers  = ((zmax2-zmin)/(Kepler*SomaSize));

% Create square packed somata within limits
xlin        = linspace(xmin,xmax2,(xmax2-xmin)/SomaSize);
ylin        = linspace(ymin,ymax2,(ymax2-ymin)/SomaSize);
zlin        = linspace(zmin,zmax2,zlayers);
[xx,yy,zz]  = meshgrid(xlin,ylin,zlin);
M           = [reshape(xx,[],1),reshape(yy,[],1),reshape(zz,[],1)];

% Shift every other z-layer to create hexagonal packing
xycells     = length(xlin)*length(ylin);
SomataGrid  = M(:,:);
for counter = 2:2:zlayers
    startpt                     = xycells*(counter-1) + 1;
    endpt                       = xycells*counter;
    SomataGrid(startpt:endpt,:) = [M(startpt:endpt,1)+(SomaSize/2),M(startpt:endpt,2),M(startpt:endpt,3)];
end

% Test whether somata lie within GCL and choose only those with centers inside GCL volume
IN_somata  = inpolyhedron(S_g.bnd,[X_g,Y_g,Z_g],SomataGrid,'FLIPNORMALS','True');
Somata     = SomataGrid(IN_somata,:);

% Save somata to file
save('Outputs/Somata_all.mat','Somata','-v7.3');
