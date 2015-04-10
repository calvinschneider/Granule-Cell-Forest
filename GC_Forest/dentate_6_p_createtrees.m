%% Select Points Within Cone Radius Distance, Rotate, and Test if Inside Elliptic Cone
% input ranges from 1 to 1186 

function dentate_6_p_createtrees(input,directory)

% Load somata and all points
load('./Outputs/Somata.mat')
load('./Outputs/Somata_pos.mat')
load('./Outputs/Points.mat')
load('./Outputs/Points_MLBinned.mat')

% Define parameters
x_split     = 50;
z_split     = 50;
subset_size = 1000;
seed        = 1;

% Define limits
xmin        = -3400;
xmax        = 3400;
zmin        = -900;
zmax        = 1700;
x_bins      = (xmax-xmin)/x_split;
z_bins      = (zmax-zmin)/z_split;

% Choose somata range based on input
input2      = str2double(input);
start_soma  = 1+(input2-1)*subset_size;
if input2*subset_size < length(Somata)
    end_soma  = input2*subset_size;
else
    end_soma  = length(Somata);
end
range           = start_soma:end_soma;

% Allocate variables
R1s             = cell(1,(end_soma-start_soma+1));
R2s             = cell(1,(end_soma-start_soma+1));
trees           = cell(1,(end_soma-start_soma+1));
counter1 = 1;

for counter = range
    % Set random seed for reproducibility
    rng(counter*seed)

    % Determine soma position
    Superficial     = Somata_pos(counter,1);
    Suprapyramidal  = Somata_pos(counter,2);
    
    % Select for ML points within 500 micron from soma
    bin_number_soma     = ceil((Somata(counter,1)-xmin)/x_split);
    bin_start_surface   = bin_number_soma - 500/x_split;
    bin_end_surface     = bin_number_soma + 500/x_split;
    if bin_start_surface < 1
        bin_start_surface = 1;
    end
    if bin_end_surface > x_bins
        bin_end_surface = x_bins;
    end
    ML_current = vertcat(ML{bin_start_surface:bin_end_surface});
    
    % Find point on ML surface closest to soma and find direction
    [k,d]           = dsearchn(ML_current,Somata(counter,:));
    closestpt       = ML_current(k,:);
    direction       = [closestpt(1,1)-Somata(counter,1),closestpt(1,2)-Somata(counter,2),closestpt(1,3)-Somata(counter,3)];
    unit_direction  = direction/norm(direction);
    
    % Define elliptical cone height and radii parameters based on soma position
    Cone_Height     = d;
    if Suprapyramidal == 1
        if Superficial == 1
            T_radius_mean   = 271;
            T_radius_stdev  = 50;
            stem_lambda     = 2.7;
            bf              = 1.35;
            NodesPerCell    = 31.5;
            NodesStDev      = 4.6;
            L_radius_mean   = 107.5;
            L_radius_stdev  = 28; 
        else
            T_radius_mean   = 182;
            T_radius_stdev  = 37;
            stem_lambda     = 0.97;
            bf              = 0.9;
            NodesPerCell    = 37.1;
            NodesStDev      = 4.8;
            L_radius_mean   = 127.5;
            L_radius_stdev  = 26;   
        end
    elseif Suprapyramidal == 0
        stem_lambda     = 1.05;    
        if Superficial == 1
            NodesPerCell    = 27.8;
            NodesStDev      = 3.7;  
            T_radius_mean   = 216;
            T_radius_stdev  = 38;        
            bf              = 1.22;   
            L_radius_mean   = 110.5;
            L_radius_stdev  = 12;   
        else
            NodesPerCell    = 27.4;
            NodesStDev      = 3.7;    
            T_radius_mean   = 143;
            T_radius_stdev  = 54;
            bf              = 1.22;
            L_radius_mean   = 101.5;
            L_radius_stdev  = 5;  
        end
    end
      
    % Choose number of stems based on truncated Poisson distribution
    while true 
        nstems = poissrnd(stem_lambda,1,1);
        if nstems >= 1 && nstems <= 4
            break
        end
    end
        
    % Choose radii from standard deviation with transverse > longitudinal
    while true 
	T_radius        = normrnd(T_radius_mean,T_radius_stdev);
	D_radius        = normrnd(L_radius_mean,L_radius_stdev);
        if T_radius > D_radius
            break
        end
    end
    
    % Create scaled elliptical cone
    vlin            = linspace(0,2*pi,100);
    ulin            = linspace(0,Cone_Height,20);
    [u,v]           = meshgrid(ulin,vlin);
    x_cone          = T_radius*sin(v).*(Cone_Height-u)/Cone_Height;
    y_cone          = D_radius*cos(v).*(Cone_Height-u)/Cone_Height;
    z_cone          = -u;
    z_cone          = z_cone + Cone_Height;
    cone            = [reshape(x_cone,[],1),reshape(y_cone,[],1),reshape(z_cone,[],1)];

    % Determine longitudinal direction
    u_soma              = Somata_pos(counter,3);
    v_soma              = Somata_pos(counter,4);  
    u1                  = u_soma - 0.01;
    u2                  = u_soma + 0.01;
    [x_pt1,y_pt1,z_pt1] = layer_eq_GCL_point(0,u1,v_soma);
    [x_pt2,y_pt2,z_pt2] = layer_eq_GCL_point(0,u2,v_soma);
    direction_l         = [x_pt2-x_pt1,y_pt2-y_pt1,z_pt2-z_pt1];
    unit_direction_l    = direction_l/norm(direction_l);
    
    % Rotate longitudinal axis and line up with previous vector
    r               = vrrotvec([unit_direction_l(:,1) unit_direction_l(:,2) 0],[0 1 0]);
    R1              = vrrotvec2mat(r);
    rotated_cone1   = cone * R1;
    r               = vrrotvec([unit_direction(:,1) unit_direction(:,2) unit_direction(:,3)],[0 0 1]);
    R2              = vrrotvec2mat(r);
    rotated_cone    = rotated_cone1 * R2;
    R1s{counter1}   = R1;
    R2s{counter1}   = R2;
    
    % Transpose cone to soma location
    x_cone          = rotated_cone(:,1) + Somata(counter,1);
    y_cone          = rotated_cone(:,2) + Somata(counter,2);
    z_cone          = rotated_cone(:,3) + Somata(counter,3);    
    
    % Only use points from bins within septotemporal radius
    xbin_start   =  ceil((min(x_cone)-xmin)/x_split);
    xbin_end     =  ceil((max(x_cone)-xmin)/x_split);
    zbin_start   =  ceil((min(z_cone)-zmin)/z_split);
    zbin_end     =  ceil((max(z_cone)-zmin)/z_split);    
    if xbin_start < 1
        xbin_start = 1;
    end
    if xbin_end > x_bins
        xbin_end = x_bins;
    end
    if zbin_start < 1
        zbin_start = 1;
    end
    if zbin_end > z_bins
        zbin_end = z_bins;
    end    
    
    points1 = cell(xbin_end - xbin_start,1);
    points2 = cell(xbin_end - xbin_start,1);
    points3 = cell(xbin_end - xbin_start,1);
    points4 = cell(xbin_end - xbin_start,1);
    points5 = cell(xbin_end - xbin_start,1);
    
    % Combine points from appropriate bins
    for bin = xbin_start:xbin_end
        points1{bin} = vertcat(pts{1}{bin}{zbin_start:zbin_end});
        points2{bin} = vertcat(pts{2}{bin}{zbin_start:zbin_end});    
        points3{bin} = vertcat(pts{3}{bin}{zbin_start:zbin_end});
        points4{bin} = vertcat(pts{4}{bin}{zbin_start:zbin_end});
        points5{bin} = vertcat(pts{5}{bin}{zbin_start:zbin_end});
    end
    
    points{1} = vertcat(points1{:});
    points{2} = vertcat(points2{:});
    points{3} = vertcat(points3{:});
    points{4} = vertcat(points4{:});    
    points{5} = vertcat(points5{:}); 
    
    % Allocate variables
    pts_subset      = cell(1,4);
    IN_econepts     = cell(1,4);
    x_pts           = cell(1,4);
    y_pts           = cell(1,4);
    z_pts           = cell(1,4);
    new_pts         = cell(1,4);
    rotated_pts     = cell(1,4);
     
    for layer = 1:5
        % Select subset of points in cubic volume      
        pts_subset{layer} = points{layer}(points{layer}(:,2)>=min(y_cone)& points{layer}(:,2)<=max(y_cone),:);
        pts_subset{layer} = pts_subset{layer}(pts_subset{layer}(:,3)>=min(z_cone)& pts_subset{layer}(:,3)<=max(z_cone),:);
        pts_subset{layer} = pts_subset{layer}(pts_subset{layer}(:,1)>=min(x_cone)& pts_subset{layer}(:,1)<=max(x_cone),:);
        
        % Untranspose points and cone
        x_pts{layer} = double(pts_subset{layer}(:,1)) - Somata(counter,1);
        y_pts{layer} = double(pts_subset{layer}(:,2)) - Somata(counter,2);
        z_pts{layer} = double(pts_subset{layer}(:,3)) - Somata(counter,3);
        new_pts{layer} = [x_pts{layer} y_pts{layer} z_pts{layer}];
        
        % Unrotate points     
        rotated_pts{layer}    = new_pts{layer} * transpose(R2);
        rotated_pts{layer}    = rotated_pts{layer} * transpose(R1);
        
        % Test if points are within elliptic cone
        rotated_pts{layer}  = [rotated_pts{layer}, (Cone_Height./rotated_pts{layer}(:,3)).^2 .* ((rotated_pts{layer}(:,1)/T_radius).^2 + (rotated_pts{layer}(:,2)/D_radius).^2)];
        IN_econepts{layer}  = pts_subset{layer}(rotated_pts{layer}(:,4)<=1 & rotated_pts{layer}(:,3)>0,:);
    end
    availablepts  = IN_econepts;
    
    % Choose percentage of nodes in each layer
    NodesN      = round(NodesPerCell + NodesStDev.*randn(1,1));
    N_GCL       = round(normrnd(0.14,0.041) * NodesN);
    N_IML       = round(normrnd(0.13,0.041) * NodesN);
    N_OML       = round(normrnd(0.115,0.00) * NodesN);    
    N_OOML      = round(normrnd(0.48,0.00) * NodesN);
    if N_OOML < 0
        N_OOML = 0;
    end
    if N_IML < 0
        N_IML = 0;
    end
    if N_GCL < 0
        N_GCL = 0;
    end
    if N_OML < 0
        N_OML = 0;
    end
    
    N_MML       = NodesN - N_OOML - N_IML - N_GCL - N_OML;

    if N_MML < 0
        N_MML = 0;
    end

    % If not enough GCL points, choose from IML
    if size(availablepts{1},1) < N_GCL
        N_IML = N_IML + (N_GCL - size(availablepts{1},1));
        N_GCL = size(availablepts{1},1);
    end

    % Choose points randomly from each layer
    pts_GCL     = availablepts{1}(randperm(size(availablepts{1}, 1), N_GCL), :);
    pts_IML     = availablepts{2}(randperm(size(availablepts{2}, 1), N_IML), :);
    pts_MML     = availablepts{3}(randperm(size(availablepts{3}, 1), N_MML), :);
    pts_OML     = availablepts{4}(randperm(size(availablepts{4}, 1), N_OML), :);    
    pts_OOML    = availablepts{5}(randperm(size(availablepts{5}, 1), N_OOML), :);  
    all_pts         = double([pts_GCL;pts_IML;pts_MML;pts_OML;pts_OOML]);
    
    % Connect chosen target points using TREES toolbox
    tree            = xMST_tree(1,[Somata(counter,1);all_pts(:,1)], ...
        [Somata(counter,2);all_pts(:,2)],...
        [Somata(counter,3);all_pts(:,3)],...
        [], [], bf, 500, [], nstems, '-b');
    % Set uniform diameter for now
    tree.D          = tree.D * 0 + 3;  
    
    trees{counter1} = tree;
    counter1        = counter1 + 1; 
end

% Save angles to file for later analysis
angles = cell(2,1);
angles{1} = R1s;
angles{2} = R2s; 

% Save files
save(sprintf('%s/Angles/%s.mat',directory,input),'angles','-v7.3')
save_tree(trees,sprintf('%s/Trees/%s.mtr',directory,input));