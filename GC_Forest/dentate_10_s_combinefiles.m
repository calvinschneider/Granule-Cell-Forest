%% Combines Parallel Files into Single File

function dentate_10_s_combinefiles(directory)

% Combine statistics files
files = dir(sprintf('%s/Stats_Tapered/*.mat',directory));
for i = 1:size(files)
    file    = sprintf('%s/Stats_Tapered/%i.mat',directory,i);
    load(file);
    f       = fieldnames(Stats);
    for j = 1:size(f,1)
        Stat.(f{j}){i} = Stats.(f{j});
    end
end
for i = 1:size(f,1)
    Stats_All.(f{i}) = vertcat(Stat.(f{i}){:});
end
save(sprintf('Outputs/Stats.mat'),'Stats_All','-v7.3')  


% Combine contraction files
files   = dir(sprintf('%s/Contraction/*.mat',directory));
M       = cell(size(files),1);
for i = 1:size(files)
    file    = sprintf('%s/Contraction/%i.mat',directory,i);
    load(file);
    M{i}    = contraction;
end
for i = 1:size(f,1)
    contraction_all = vertcat(M{:});
end
save(sprintf('Outputs/Contraction.mat'),'contraction_all','-v7.3')  


% Combine diameter files
files   = dir(sprintf('%s/Diameters/*.mat',directory));
M       = cell(1186,1);
for i = 1:size(files)
    file    = sprintf('%s/Diameters/%i.mat',directory,i);
    load(file);
    M{i}    = diameters;
end
diameters_all = vertcat(M{:});

% Save diameter information based on 10 micron euclidean distance bins
all = vertcat(diameters_all{:});
length(all)
bin_size = 10;
[bincounts,index] = histc(all(:,1),0:bin_size:max(all(:,1)));
means = zeros(length(bincounts),3);
for i = 1:length(bincounts)
    means(i,1) = (i-1) * bin_size;
    means(i,2) = mean(all(index == i,2));
    means(i,3) = std(all(index == i,2));
end
save(sprintf('Outputs/Diameters.mat'),'means','-v7.3')

