This code was used in the publication:
Schneider CJ, Cuntz H, Soltesz I
Linking Macroscopic with Microscopic Neuroanatomy using Synthetic Neuronal Populations
PLOS Comp Biol

The code was run on a unix cluster with the SGE Queueing System. 
The generation process is broken down into steps in order to be modified easily.
Scripts to submit each step to SGE can be found in /Submit for each step.
Each step in the generation process was first compiled into executables with the MATLAB mcc compiler (see compile.sh)
Parallel files are saved in their own directory (which should be changed by the particular user) and then later combined into single files in the /Outputs folder.
"s" - denotes a step run on a single processor (serial)
"p" - denotes a step run using parallel computing
TREES_Analysis contains files from/modified from the TREES Toolbox http://www.treestoolbox.org/
Supporting_Files contains other files necessary to run the generation process


The output directory for files created in parallel and the queues available should be changed in all scripts in the Submit/ folder before use.

Generation Process:
dentate_1_s_fillsomata 		- Creates a grid of hexagonal packed somata and outputs somata with centers inside of GCL
dentate_2_p_somataposition	- Determines whether somata are infra/suprapyramidal and deep/superficial and eliminates those within a cell radius of the GCL border
dentate_3_s_combinesomata	- Combines somata and position files into single files for later retrieval
dentate_4_p_findallpoints	- Finds all possible target points in the GCL and ML at 1 micron interval (must be run 5 times total, for sublayers 1-5 in submit script)
dentate_5_s_choosepoints 	- Chooses a subset of points based on input laminar ratio
dentate_6_p_createtrees 	- Chooses target points and creates trees
dentate_7_p_jitter          	- Takes created trees and adds spatial jitter
dentate_8_p_taper           	- Takes jittered trees and adds diameter taper
dentate_9_p_analysis		- Analyzes generated GC population properties
dentate_10_s_combinefiles	- Combines analysis files created in parallel into single files
dentate_11_s_printstats     	- Prints out Tables 2 and 3 from paper (run interactively, not through SGE)

Outputs:
Contraction.mat             - Contraction values for all branches (used to create histogram in Figure 3)
Diameters.mat               - Binned averages for diameters at each 10 micron interval euclidean distance from soma (used for plot in Figure 3)
Points.mat                  - All available points (used in generation)
Points_MLBinned.mat         - Outer molecular layer boundary points (used in generation)
Somata_all.mat              - All somata with centers inside of GCL (before trimming those that are not fully contained in GCL)
Somata_pos.mat              - Position information for all somata used
Somata.mat                  - 3D coordinates of somata used
Stats.mat                   - Statistics for all trees

All generated trees are located in your specified directory, in the Trees_Tapered folder. 