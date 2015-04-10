% Analyze Stats Based on Position in Granule Cell Layer

function dentate_9_p_analysis(input,directory)

global trees
trees       = {};
subset_size = 1000;

% Load Somata
load('Outputs/Somata.mat')
load('Outputs/Somata_pos.mat')


% Load trees and angles associated with input
load(sprintf('%s/Angles/%s.mat',directory,input))
load_tree(sprintf('%s/Trees_Tapered/%s.mtr',directory,input))
trees       = trees{1};
input2      = str2double(input);
start_soma  = 1+(input2-1)*subset_size;
if input2*subset_size < length(Somata)
    end_soma  = input2*subset_size;
else
    end_soma  = length(Somata);
end
range           = start_soma:end_soma;    

% Preallocate variables
Stats = struct;

% Get stats from TREES toolbox
s = stats_tree(trees,[],[],'-x');

% Rotate trees and get spanning information
trees_rotated   = rotate_created(trees,angles{1},angles{2}); 
[spanning,~]    = gscale_tree2(trees_rotated,' ');

% Calculate dendritic segments for each tree and append to matrix
% Allocate variables
Stats.type          = uint8(Somata_pos(range,1:2));
Stats.n_dend        = zeros(length(trees),1,'uint8');
Stats.n_branches    = zeros(length(trees),1,'uint8');
Stats.max_BO        = zeros(length(trees),1,'uint8');
Stats.BO            = cell(length(trees),1);
Stats.tip_dist      = zeros(length(trees),1);
Stats.tip_path      = zeros(length(trees),1);
Stats.sholl         = cell(length(trees),1);
Stats.asym          = zeros(length(trees),1);
Stats.inter_plen    = zeros(length(trees),1);
Stats.term_plen     = zeros(length(trees),1);
Stats.percents      = zeros(length(trees),16,'single');
Stats.length        = single(s.gstats.len);
Stats.spread_trans  = single(spanning.xlims{1});
Stats.height        = single(spanning.ylims{1});
Stats.spread_long   = single(spanning.zlims{1});
BO_all              = cell(length(trees),1);
B                   = cell(4,1);
T                   = cell(4,1);
X                   = cell(4,1);
Y                   = cell(4,1);
Z                   = cell(4,1);
S1                  = cell(4,1);
all_B               = cell(4,1);
all_T               = cell(4,1);
B_pts               = cell(length(trees),1);
T_pts               = cell(length(trees),1);
B{1}                = cell(length(trees),1);
B{2}                = cell(length(trees),1);
B{3}                = cell(length(trees),1);
B{4}                = cell(length(trees),1);
T{1}                = cell(length(trees),1);
T{2}                = cell(length(trees),1);
T{3}                = cell(length(trees),1);
T{4}                = cell(length(trees),1);
eucdist             = cell(length(trees),1);
pathdist            = cell(length(trees),1);
contraction         = cell(length(trees),1);
branch_indices      = cell(length(trees),1);

for sublayer = 1:4
    % Define surface for GCL and ML thirds
    if sublayer == 1
        [x_i,y_i,z_i]       = layer_eq_GCL(-1.95);
        [x_m,y_m,z_m]       = layer_eq_GCL(-1.0);
        [x_o,y_o,z_o]       = layer_eq_GCL(0);
        X{sublayer}         = [x_o;x_m;x_i];
        Y{sublayer}         = [y_o;y_m;y_i];
        Z{sublayer}         = [z_o;z_m;z_i];
        [~,S1{sublayer}]    = alphavol([X{sublayer}(:),Y{sublayer}(:),Z{sublayer}(:)],120);    
    elseif sublayer > 1
        [x_i,y_i,z_i]       = layer_eq_ML(sublayer-2);
        [x_o,y_o,z_o]       = layer_eq_ML(sublayer-1);
        X{sublayer}         = [x_o;x_i];
        Y{sublayer}         = [y_o;y_i];
        Z{sublayer}         = [z_o;z_i];
        [~,S1{sublayer}]    = alphavol([X{sublayer}(:),Y{sublayer}(:),Z{sublayer}(:)],150);
    end
end

% Find statistics for individual cells
for i = 1:length(trees)
    % Define topological points and remove soma branch points
    Soma        = Somata(range(i),:);
    iB          = B_tree(trees{i});
    iT          = T_tree(trees{i}); 
    B_pts{i}    = [trees{i}.X(iB) trees{i}.Y(iB) trees{i}.Z(iB)];
    T_pts{i}    = [trees{i}.X(iT) trees{i}.Y(iT) trees{i}.Z(iT)];

    % Find asymmetry
    asym                = asym_tree(trees{i});
    Stats.asym(i)       = mean(asym(iB));
    
    % Remove soma as a branch point
    if ismember(Soma,B_pts{i},'rows') > 0
        iB(1,1)     = 0;
        B_pts{i}(1,:)  = [];
    end
    iBT     = iB | iT;
    
    Stats.n_dend(i,1)     = size(find(PL_tree(trees{i})==1),1);
    Stats.n_branches(i,1) = size(B_pts{i},1)*2 + Stats.n_dend(i,1);
               
    % Branch order
    BO_all{i}           = BO_tree(trees{i});
    Stats.BO{i}         = uint8(BO_all{i}(iBT));
    Stats.max_BO(i,1)   = max(Stats.BO{i})+1;
    
    % Find euclidean and path distance to tips
    eucl                = eucl_tree(trees{i});
    Plen                = Pvec_tree(trees{i});
    Stats.tip_dist(i)   = mean(eucl(iT));
    Stats.tip_path(i)   = mean(Plen(iT));
    
    % Find sholl distribution
    diam                = 50;
    Stats.sholl{i}      = uint8(sholl_tree(trees{i},diam));
    
    % Break down tree into branches
    [sect,~]        = dissect_tree(trees{i});
    len             = len_tree(trees{i});
    sect            = sect(sect(:,1) ~= sect(:,2),:);
    n_branches      = length(sect);
    % Get indices for each branch
    ipar            = ipar_tree(trees{i});    
    
    % Allocate variables
    eucdist{i}           = zeros(n_branches,1);
    pathdist{i}          = zeros(n_branches,1);
    contraction{i}       = zeros(n_branches,1,'single');
    branch_indices{i}    = zeros(n_branches,1);
    
    for j = 1:n_branches
        % Determine all nodes for each branch from ipar
        start               = sect(j,1);
        stop                = sect(j,2);
        
        parent_line         = ipar(stop,:);
        stop_index          = find(parent_line == start);
        if isempty(stop_index) == 1
            stop_index = find(parent_line > 0,1,'last')+1;
        end
        branch_indices{j}   = transpose(parent_line(1:stop_index-1));
        
        % Get euclidean distance from branch beginning and end
        startxyz                = [trees{i}.X(start) trees{i}.Y(start) trees{i}.Z(start)];    
        stopxyz                 = [trees{i}.X(stop) trees{i}.Y(stop) trees{i}.Z(stop)];    
        eucdist{i}(j,1)         = pdist([startxyz;stopxyz],'euclidean');
        pathdist{i}(j,1)        = sum(len(branch_indices{j}));
        contraction{i}(j,1)     = eucdist{i}(j,1)/pathdist{i}(j,1);
    end    
    
    % Get intermediate and terminal branch lengths
    sect                = dissect_tree(trees{i});
    blen                = diff (Plen(sect), [], 2);
    [~,term_i]          = intersect(sect(:,2),find(iT));
    [~,inter_i]         = setdiff(sect(:,2),find(iT));
    inter_plen1         = blen(inter_i);
    Stats.inter_plen(i) = mean(inter_plen1(inter_plen1>0.01));
    Stats.term_plen(i)  = mean(blen(term_i));
end

all_B_pts = vertcat(B_pts{:});
all_T_pts = vertcat(T_pts{:});

% Determine which layer branch and termination points are in
for sublayer = 1:4
    IN_B            = inpolyhedron(S1{sublayer}.bnd,[X{sublayer},Y{sublayer},Z{sublayer}],all_B_pts,'FLIPNORMALS','True');
    IN_T            = inpolyhedron(S1{sublayer}.bnd,[X{sublayer},Y{sublayer},Z{sublayer}],all_T_pts,'FLIPNORMALS','True');
    all_B{sublayer} = all_B_pts(IN_B,:);
    all_T{sublayer} = all_T_pts(IN_T,:);
end

% Write out branch and termination laminar distribution percentages
for i = 1:length(trees)
    for sublayer = 1:4
       B{sublayer}{i} = intersect(B_pts{i},all_B{sublayer},'rows');
       T{sublayer}{i} = intersect(T_pts{i},all_T{sublayer},'rows');
    end
    B_Total = size(B{1}{i},1) + size(B{2}{i},1) + size(B{3}{i},1) + size(B{4}{i},1);
    T_Total = size(T{1}{i},1) + size(T{2}{i},1) + size(T{3}{i},1) + size(T{4}{i},1);

    % Write out percentages to matrix
    Stats.percents(i,:) = [100*size(B{1}{i},1)/B_Total 100*size(B{2}{i},1)/B_Total 100*size(B{3}{i},1)/B_Total 100*size(B{4}{i},1)/B_Total...
                     100*size(T{1}{i},1)/B_Total 100*size(T{2}{i},1)/T_Total 100*size(T{3}{i},1)/T_Total 100*size(T{4}{i},1)/T_Total...
                     100*(size(B{1}{i},1)+size(T{1}{i},1))/(B_Total+T_Total) 100*(size(B{2}{i},1)+size(T{2}{i},1))/(B_Total+T_Total) 100*(size(B{3}{i},1)+size(T{3}{i},1))/(B_Total+T_Total) 100*(size(B{4}{i},1)+size(T{4}{i},1))/(B_Total+T_Total)...
                     size(B{1}{i},1)+size(T{1}{i},1) size(B{2}{i},1)+size(T{2}{i},1) size(B{3}{i},1)+size(T{3}{i},1) size(B{4}{i},1)+size(T{4}{i},1)];
end

% Save outputs
save(sprintf('%s/Stats_Tapered/%s.mat',directory,input),'Stats','-v7.3')
save(sprintf('%s/Contraction/%s.mat',directory,input),'contraction','-v7.3')
