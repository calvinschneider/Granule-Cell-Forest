%% Combines Parallel Soma Files

function dentate_3_s_combinesomata(directory)

% Combine Soma Position Files
% Determine number of files
files       = dir(sprintf('%s/Soma_Position/*.mat',directory));

% Concatenate positions of all somata and write to file
pos = cell(length(files),1);
for i = 1:length(files)
    file        = sprintf('%s/Soma_Position/%i.mat',directory,i);
    load(file);
    pos{i}      = Somata_pos;
end
Somata_pos = vertcat(pos{1:length(files)});
save('Outputs/Somata_pos.mat','Somata_pos','-v7.3');


% Combine Trimmed Soma Files
% Cycle through folder and determine number of files
files = dir(sprintf('%s/Soma_Trimmed/*.mat',directory));

% Concatenate all somata and write to file
trimmed = cell(length(files),1);
for i = 1:length(files)
    file        = sprintf('%s/Soma_Trimmed/%i.mat',directory,i);
    load(file);
    trimmed{i}  = Somata;
end
Somata = vertcat(trimmed{1:length(files)});
save('Outputs/Somata.mat','Somata','-v7.3');
