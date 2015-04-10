%% Adds Taper to Generated Trees

function dentate_8_p_taper(input,directory)

% Define variables
global trees
trees       = {};
subset_size = 1000;
seed        = 1;

% Load Files
load('Outputs/Somata.mat')
load('TREES1.15/construct/P.mat')
load('TREES1.15/construct/ldend.mat')

% Choose tree range based on input
load_tree(sprintf('%s/Trees_Jittered/%s.mtr',directory,input))
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
diameters   = cell(length(trees),1);
counter1 = 1;

for counter = range
    % Set random number generator for reproducible results
    rng(counter*seed)
    
    % Choose taper parameters and make non-negative
    while true
        scale           = normrnd(0.0275,0.020);
        offset          = normrnd(0.66,0.105);
        root_distance   = normrnd(30,5);
        if scale > 0 && offset > 0 && root_distance > 0
            break
        end
    end
    trees2{counter1} = quadfuncdiam_tree(trees{counter1},[scale,0.036,offset,root_distance],@(x) exp(x)-1,'',P,ldend);
    
    % Write out euclidean distance and diameter for later analysis
    diameters{counter1}        = single(zeros(size(trees2{counter1}.D,1),2));
    diameters{counter1}(:,1)   = eucl_tree(trees2{counter1});
    diameters{counter1}(:,2)   = trees2{counter1}.D;
    counter1 = counter1 + 1;
end
trees = trees2;

% Save outputs
save_tree(trees,sprintf('%s/Trees_Tapered/%s.mtr',directory,input));    
save(sprintf('%s/Diameters/%s.mat',directory,input),'diameters','-v7.3')   
