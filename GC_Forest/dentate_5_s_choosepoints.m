%% Select Random Points within the GCL and Molecular Layer and Keep in Bins

function dentate_5_s_choosepoints(directory)

% Define subset of points and laminar distribution
TotalPoints  = 200000000;
Nodes_GCL   = round(0.14 * TotalPoints);
Nodes_IML   = round(0.14 * TotalPoints);
Nodes_MML   = round(0.12 * TotalPoints);
Nodes_OML   = round(0.20 * TotalPoints);
Nodes_OOML  = round(0.40 * TotalPoints);

% Define parameters
zmin        = -900;
zmax        = 1700;
z_split     = 50;
z_bins      = (zmax - zmin)/z_split;
seed        = 1;

% Define thirds
layers = {'GCL','IML','MML','OML','OOML'};
pts = cell(1,length(layers));

for k = 1:length(layers)
    
    % Set random number generator for reproducible results
    rng(k*seed)
    
    % Cycles through available .mat files and determines number of points
    n               = eval(['Nodes_' sprintf('%s',layers{k})]);
    files           = dir(sprintf('%s/All_Points/%s/*.mat',directory,layers{k}));
    points_per_bin  = cell(size(files,1),1);
    sum_bins        = cell(size(files,1),1);
    all_pts         = cell(size(files,1),1);
    bin_totals      = zeros(size(files,1),1);
    sum             = 0;
    for i = 1:size(files,1)
        load(sprintf('%s/All_Points/%s/%i.mat',directory,layers{k},i));
        all_pts{i}  = M;
        for j = 1:z_bins
            points_per_bin{i}(j,1) = size(all_pts{i}{j},1);
            sum                    = sum + size(all_pts{i}{j},1);
            sum_bins{i}(j,1)       = sum;
        end
        bin_totals(i,1) = sum;
    end

    % Choose random numbers in correct range
    totalpoints     = sum;
    point_numbers   = int64(reshape(randperm(totalpoints,n),[],1));
    
    % Select points from appropriate files
    pts{k}          = cell(size(all_pts,1),1);
    for i = 1:size(all_pts,1)
        pts{k}{i}   = cell(z_bins,1);
        if i == 1
            point_numbers1 = point_numbers(point_numbers > 0 & point_numbers < bin_totals(i,1),1);
        else
            point_numbers1 = point_numbers(point_numbers > bin_totals(i-1,1) & point_numbers < bin_totals(i,1),1);
        end
        
        % Extract chosen corresponding points from loaded points
        for j = 1:z_bins
            if j == 1
                if i == 1
                    if sum_bins{i}(j,1) > 0
                        points_to_use   = int64(point_numbers1(point_numbers1(:,1)> 0 & point_numbers1(:,1) < sum_bins{i}(j,1),1));
                        pts{k}{i}{j}    = int16(all_pts{i}{j}(points_to_use,:));
                    end
                else
                    if sum_bins{i}(j,1) > sum_bins{i-1}(z_bins,1)
                        points_to_use   = int64(point_numbers1(point_numbers1(:,1)> sum_bins{i-1}(z_bins,1) & point_numbers1(:,1) < sum_bins{i}(j,1),1) - sum_bins{i-1}(z_bins,1));
                        pts{k}{i}{j}    = int16(all_pts{i}{j}(points_to_use,:));
                    end
                end
            else
                if sum_bins{i}(j,1) > sum_bins{i}(j-1,1)
                    points_to_use   = int64(point_numbers1(point_numbers1(:,1)> sum_bins{i}(j-1,1) & point_numbers1(:,1) < sum_bins{i}(j,1),1) - sum_bins{i}(j-1,1));
                    pts{k}{i}{j}    = int16(all_pts{i}{j}(points_to_use,:));
                end
            end
        end
    end
        
end

% Define outer ML surface for later orientation of cones
[x_ml,y_ml,z_ml]    = layer_eq_ML_2(3);
ML_surface          = [x_ml,y_ml,z_ml];

% Define ML parameters
split_x     = 50;
xmin        = floor(min(x_ml)/split_x)*split_x;
xmax        = ceil(max(x_ml)/split_x)*split_x;
x_bins      = (xmax-xmin)/split_x;
ML          = cell(1,x_bins);

% Split outer ML surface into septotemporal bins for speed
for bin_x = 1:x_bins
    ML{bin_x}         = ML_surface(ML_surface(:,1)>=(xmin+(bin_x-1)*split_x+1)& ML_surface(:,1)<=(xmin+bin_x*split_x),:);
end

save('Outputs/Points.mat','pts','-v7.3');
save('Outputs/Points_MLBinned.mat','ML','-v7.3')
