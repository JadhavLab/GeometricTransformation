+-------------------------------------------------------------------------+
|      Geometric transformation of cognitive maps for generalization      |
|                    across hippocampal-prefrontal circuits               |
+-------------------------------------------------------------------------+
README.txt
Copyright (C) 2023, Wenbo Tang, version 1.0
All rights reserved.



BRIEF
=====

Code accompanying the paper: Tang, W., Shin, J. D. & Jadhav, S. P. (2023). Geometric transformation of cognitive maps for generalization across hippocampal-prefrontal circuits. Cell Reports.


GETTING STARTED
===============

Launch MATLAB and cd into the directory containing the code (e.g. '/CellRep_2023/').


Time-filter framework: Jadhav Lab (Brandeis University) and Frank Lab (UCSF)
Other files in the directory (with all sub-folders) needed in path for running the time-filter framework and other analysis:
\usrlocal\   
\Src_Matlab\  
\Sleep_Code\

Toolboxes required:
- Libsvm (version 3.12; https://www.csie.ntu.edu.tw/~cjlin/libsvm/) 

- Uniform Manifold Approximation and Projection (UMAP) (version 4.1; https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)

- MatPlotLib (version 2.1.3; https://www.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps)


These codes were originally created in the MATLAB 2017a. All main scripts have plotting functions built-in to generate the figures shown in the paper.


FILES and FOLDERS
=================
  ./Figure1
  cal_behavperform.m	: main script for calculating performance for each animal
  behavperform_gather.m	: main script for generating statistics of performance in Figure 1
  
  ./Figure1/Performance_files
  subfolder containing all performance files generated by cal_behavperform.m


  ./Figure2
  Plot_linearRateMaps_allcells.m  : script for plotting the sorted linearized rate maps in  Figure 2 
  cal_PVsimilarity.m	: compute PV similarity across environments
  decoding_position_novelfamiliar_batch.m  : main script for decoding animal's current  position using rate maps
  decoding_position_novelfamiliar_templateN2_batch.m  : main script for decoding animal's current position using rate maps from N’
  decoding_position_novelfamiliar_trial_batch.m  : main script for get decoding error of animal's current position using rate maps from N' in a trial-by-trial basis
 cal_confusionMat_novelfamiliar.m    : script for plotting the confusion matrix in Figure  2

 ./Figure2/Decodepos
  subfolder containing all decoding results generated by decoding_position_novelfamiliar_batch.m and decoding_position_novelfamiliar_templateN2_batch.m


  ./Figure3
  cal_SI.m		: main script for calculating trajectory selectivity index
  Plot_linfields_sortedbySI.m  : script for plotting plotting all linearized rate maps sorted by trajectory selectivity index in  Figure 3  
  cal_UMAP.m		: main script for UMAP transformation of neuronal population activity
  cal_UMAP_TrajPhase_vs_Spatial_distance.m : main script for calculating the distance of neural states of the same spatial location vs. the same trajectory phase on UMAP neural manifolds
  cal_UMAP_INSeq_distance.m : main script for calculating INSeq vs. OUTSeq trajectory distance based on UMAP manifolds of neural population activity
  cal_UMAP_FNdistance_shuffle.m : main script for calculating the distance of neural states between N' and shuffled neural manifolds
  cal_UMAP_FNdistance.m : main script for calculating the distance of neural states between N' and F neural manifolds, and comparing to the shuffles

  ./Figure3/Supplemental : 
  subfolder containing supplemental analysis related to Figure 3
  cal_path_equivalence.m : main script for calculating path-equivalent coefficient
  cal_PV_TrajPhase_vs_Spatial_distance.m : main script for calculating the distance of neural states of the same spatial location vs. the same trajectory phase in the original state space
  cal_PV_INSeq_distance.m : main script for calculating INSeq vs. OUTSeq trajectory distance in the original state space
  cal_PV_FNdistance_shuffle.m : main script for calculating the distance of neural states between N' and shuffled neural activity in the original state space
  cal_PV_FNdistance.m : main script for calculating the distance of neural states between N' and F neural activity in the original state space, and comparing to the shuffles
  
  ./Figure3/SingleTrial_ratemaps : 
  subfolder containing all files for single-trial firing rates


  ./Figure4
  cal_all_dichotomies.m	: main script for getting all dichotomies (clusters) in Figure 4
  CCGP_trajPhase.m	: main script for CCGP of trajectory phases using linear SVMs
  CCGP_taskSeq.m	: main script for CCGP of task sequences using linear SVMs
  CCGP_environment.m	: main script for CCGP of different environments using linear SVMs
  decoding_CCGP_linearSVM_simple.m: helper function that computers CCGP using linear SVMs, and tests significance using trial-label shuffles 

  ./Figure4/Supplemental : 
  subfolder containing supplemental analysis related to Figure 4
  decoding_dichotomy_main.m	: main script for decoding all dichotomies using 4-fold cross-validation
  decoding_dichotomy_linearSVM.m: helper function that computers decoding accuracy using linear SVMs, and tests significance using trial-label shuffles 


CITING OUR WORK
===============

If you find the code useful, please cite the code source and the paper:
    Tang, W., Shin, J. D. & Jadhav, S. P. (2023). Geometric transformation of cognitive maps for generalization across hippocampal-prefrontal circuits. Cell Reports.


CONTACT
=======
Bug reports, comments and questions are appreciated.
Please write to: 
	Wenbo Tang <wenbo.tang07@gmail.com>
