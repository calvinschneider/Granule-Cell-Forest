%% Finds All Available Points Inside GCL and ML
% sublayer ranges from 1 to 5, input ranges from 1 to 136 for each sublayer

function dentate_4_p_findallpoints(sublayer,input,directory)

% Choose 50 micron bin in x direction based on input number
x       = str2double(input);
x_split = 50;

% Split z direction into 50 micron bins for speed when loading in choosing target point step
z_split = 50;

% Change folder name based on sublayer input
sublayer2       = str2double(sublayer);
if sublayer2 == 1
    layer_name = 'GCL';
elseif sublayer2 == 2
    layer_name = 'IML';
elseif sublayer2 == 3
    layer_name = 'MML';
elseif sublayer2 == 4
    layer_name = 'OML';
elseif sublayer2 == 5
    layer_name = 'OOML';
end

% Define volume to test points
xmin        = -3400;
xmax        = 3400;
ymin        = -800;
ymax        = 4500;
zmin        = -900;
zmax        = 1700;
z_bins      = (zmax - zmin)/z_split;

% Define surfaces
if sublayer2 == 1
    % GCL
    [x_i,y_i,z_i]   = layer_eq_GCL(-1.95);
    [x_m,y_m,z_m]   = layer_eq_GCL(-1);
    [x_o,y_o,z_o]   = layer_eq_GCL(0);
    X               = [x_o;x_m;x_i];
    Y               = [y_o;y_m;y_i];
    Z               = [z_o;z_m;z_i];
    [~,S1]          = alphavol([X(:),Y(:),Z(:)],120);   
elseif sublayer2 > 1 && sublayer2 < 4
    % IML and MML
    [x_i,y_i,z_i]   = layer_eq_ML(sublayer2-2);
    [x_o,y_o,z_o]   = layer_eq_ML(sublayer2-1);
    X               = [x_o;x_i];
    Y               = [y_o;y_i];
    Z               = [z_o;z_i];
    [~,S1]          = alphavol([X(:),Y(:),Z(:)],150);
elseif sublayer2 == 4
    % Inner 3/4 of OML
    [x_i,y_i,z_i]   = layer_eq_ML(2);
    [x_o,y_o,z_o]   = layer_eq_ML(2.75);
    X               = [x_o;x_i];
    Y               = [y_o;y_i];
    Z               = [z_o;z_i];
    [~,S1]          = alphavol([X(:),Y(:),Z(:)],150);
elseif sublayer2 == 5
    % Outer 1/4 of OML
    [x_i,y_i,z_i]   = layer_eq_ML(2.75);
    [x_o,y_o,z_o]   = layer_eq_ML(3);
    X               = [x_o;x_i];
    Y               = [y_o;y_i];
    Z               = [z_o;z_i];
    [~,S1]          = alphavol([X(:),Y(:),Z(:)],150);
end

% Create evenly spaced square grid of points 1 micron apart and test if in surface 
xlin            = linspace(xmin+1+(x-1)*x_split,xmin+x*x_split,x_split);
M               = cell(z_bins,1);

for z = 1:z_bins
    % Test square grid of points
    ylin        = linspace(ymin+1,ymax,(ymax-ymin));
    zlin        = linspace(zmin+1+(z-1)*z_split,zmin+(z*z_split),z_split);
    [xx,yy,zz]  = meshgrid(xlin,ylin,zlin);
    test_pts    = [reshape(xx,[],1),reshape(yy,[],1),reshape(zz,[],1)];
    IN          = inpolyhedron(S1.bnd,[X,Y,Z],test_pts,'FLIPNORMALS','True');
    M{z}        = int16(test_pts(IN,:));
end
save(sprintf('%s/All_Points/%s/%s.mat',directory,layer_name,input),'M','-v7.3')
