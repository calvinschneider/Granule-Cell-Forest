module load MATLAB/r2014a

mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_1_s_fillsomata.m      -a ./Supporting_Files
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_2_p_somataposition.m  -a ./Supporting_Files
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_3_s_combinesomata.m
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_4_p_findallpoints.m   -a ./Supporting_Files
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_5_s_choosepoints.m    -a ./Supporting_Files
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_6_p_createtrees.m     -a ./Supporting_Files  -a ./TREES1.15
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_7_p_jitter.m          -a ./Supporting_Files  -a ./TREES1.15
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_8_p_taper.m           -a ./Supporting_Files  -a ./TREES1.15 
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_9_p_analysis.m        -a ./Supporting_Files  -a ./TREES1.15
mcc -m -R -singleCompThread -R '-nodisplay,-nojvm' dentate_10_s_combinefiles.m
