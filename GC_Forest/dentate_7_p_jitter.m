%% Adds Jitter to Generated Trees
% input ranges from 1 to 1186

function dentate_7_p_jitter(input,directory)

% Define variables
global trees
trees       = {};
subset_size = 1000;
seed        = 1;

% Load somata
load('Outputs/Somata.mat')
load('Outputs/Somata_pos.mat')

% Choose tree range based on input
load_tree(sprintf('%s/Trees/%s.mtr',directory,input))
trees = trees{1};
input2      = str2double(input);
start_soma  = 1+(input2-1)*subset_size;
if input2*subset_size < length(Somata)
    end_soma  = input2*subset_size;
else
    end_soma  = length(Somata);
end
range       = start_soma:end_soma;    

% Allocate variables
trees2      = cell(length(trees),1);
counter1    = 1;
for counter = range
    % Set random number generator for reproducible results
    rng(counter*seed)
    
    % Resample tree at 5 micron interval
    tree = resample_tree(trees{counter1},5,' ');

    % Define spatial jitter amplitude
    stde_mean   = 0.4;
    stde_stdev  = 0.1;
    stde_mean2   = 0.25;
    stde_stdev2  = 0.125;
    
    % Jitter tree
    stde  = normrnd(stde_mean,stde_stdev);
    stde2 = normrnd(stde_mean2,stde_stdev2);
    tree2 = jitter_tree2(tree,stde,10,' ');
    trees2{counter1} = jitter_tree2(tree2,stde2,2,' ');
    counter1 = counter1 + 1;
end
trees = trees2;

% Save jittered trees
save_tree(trees,sprintf('%s/Trees_Jittered/%s.mtr',directory,input));
