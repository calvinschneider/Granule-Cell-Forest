%% Recreates Tables 2 and 3 From Paper

load('Outputs/Stats.mat')

% Split trees into subgroups
n_groups        = 7;
subgroup        = cell(n_groups,1);
all_stats       = cell(n_groups,1);
subgroup{1}     = find(Stats_All.type(:,1)>=0);                             % all trees
subgroup{2}     = find(Stats_All.type(:,2)==1);                             % suprapyramidal
subgroup{3}     = find(Stats_All.type(:,2)==0);                             % infrapyramidal
subgroup{4}     = find(Stats_All.type(:,1)==1 & Stats_All.type(:,2)==1);    % suprapyramidal superficial
subgroup{5}     = find(Stats_All.type(:,1)==0 & Stats_All.type(:,2)==1);    % suprapyramidal deep
subgroup{6}     = find(Stats_All.type(:,1)==1 & Stats_All.type(:,2)==0);    % infrapyramidal superficial
subgroup{7}     = find(Stats_All.type(:,1)==0 & Stats_All.type(:,2)==0);    % infrapyramidal deep 

for i = 1:n_groups
    all_stats{i} = [...
        mean(Stats_All.n_dend(subgroup{i})),     std(single(Stats_All.n_dend(subgroup{i}))),...
        mean(Stats_All.n_branches(subgroup{i})), std(single(Stats_All.n_branches(subgroup{i}))),...
        mean(Stats_All.max_BO(subgroup{i})),     std(single(Stats_All.max_BO(subgroup{i}))),...
        mean(Stats_All.spread_trans(subgroup{i},2) - Stats_All.spread_trans(subgroup{i},1)),...
        std( Stats_All.spread_trans(subgroup{i},2) - Stats_All.spread_trans(subgroup{i},1)),...
        mean(Stats_All.spread_long(subgroup{i},2)  - Stats_All.spread_long(subgroup{i},1)),...
        std( Stats_All.spread_long(subgroup{i},2)  - Stats_All.spread_long(subgroup{i},1)),...
        mean(Stats_All.length(subgroup{i})),     std(Stats_All.length(subgroup{i})),...
        mean(Stats_All.tip_path(subgroup{i})),   std(Stats_All.tip_path(subgroup{i})),...
        mean(Stats_All.inter_plen(subgroup{i})), std(Stats_All.inter_plen(subgroup{i})),...
        mean(Stats_All.term_plen(subgroup{i})),  std(Stats_All.term_plen(subgroup{i})),...
        mean(Stats_All.asym(subgroup{i})),       std(Stats_All.asym(subgroup{i}))];
end

% Table 2 Statistics
display(sprintf('# Dendrites                            %.1f +- %.1f',  all_stats{1}(1), all_stats{1}(2)));
display(sprintf('# Dendritic Branches                   %.0f +- %.0f',  all_stats{1}(3), all_stats{1}(4)));
display(sprintf('Max Branch Order                       %.1f +- %.1f',  all_stats{1}(5), all_stats{1}(6)));
display(sprintf('Transverse Spread (?m)                 %.0f +- %.0f',  all_stats{1}(7), all_stats{1}(8)));
display(sprintf('Longitudinal Spread (?m)               %.0f +- %.0f',  all_stats{1}(9), all_stats{1}(10)));
display(sprintf('Total Dendritic Length (?m)            %.0f +- %.0f',  all_stats{1}(11),all_stats{1}(12)));
display(sprintf('Mean Pathlength to Terminal Tips (?m)  %.0f +- %.0f',  all_stats{1}(13),all_stats{1}(14)));
display(sprintf('Mean Intermediate Branch Length (?m)   %.0f +- %.0f',  all_stats{1}(15),all_stats{1}(16)));
display(sprintf('Mean Terminal Branch Length (?m)       %.0f +- %.0f',  all_stats{1}(17),all_stats{1}(18)));
display(sprintf('Topological Asymmetry                  %.2f +- %.2f\n',all_stats{1}(19),all_stats{1}(20)));

% Table 3 Statistics
display(sprintf('# Dendritic Branches          Suprapyramidal               %.0f +- %.0f',all_stats{2}(3), all_stats{2}(4)));
display(sprintf('                              Infrapyramidal               %.0f +- %.0f',all_stats{3}(3), all_stats{3}(4)));
display(sprintf('Total Dendritic Length (?m)   Suprapyramidal               %.0f +- %.0f',all_stats{2}(11),all_stats{2}(12)));
display(sprintf('                              Infrapyramidal               %.0f +- %.0f',all_stats{3}(11),all_stats{3}(12)));
display(sprintf('# Dendrites                   Suprapyramidal Superficial   %.1f +- %.1f',all_stats{4}(1), all_stats{4}(2)));
display(sprintf('                              Suprapyramidal Deep          %.1f +- %.1f',all_stats{5}(1), all_stats{5}(2)));
display(sprintf('Max Branch Order              Suprapyramidal Superficial   %.1f +- %.1f',all_stats{4}(5), all_stats{4}(6)));
display(sprintf('                              Suprapyramidal Deep          %.1f +- %.1f',all_stats{5}(5), all_stats{5}(6)));
display(sprintf('Transverse Spread (?m)        Suprapyramidal Superficial   %.0f +- %.0f',all_stats{4}(7), all_stats{4}(8)));
display(sprintf('                              Suprapyramidal Deep          %.0f +- %.0f',all_stats{5}(7), all_stats{5}(8)));
display(sprintf('                              Infrapyramidal Superficial   %.0f +- %.0f',all_stats{6}(7), all_stats{6}(8)));
display(sprintf('                              Infrapyramidal Deep          %.0f +- %.0f',all_stats{7}(7), all_stats{7}(8)));