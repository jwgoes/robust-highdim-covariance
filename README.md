# robust-highdim-covariance
Code for "Robust Sparse Covariance Estimation by Thresholding Tyler’s M-Estimator" (https://arxiv.org/pdf/1706.08020.pdf).


This file describes three matlab functions that are used in our paper and how to use them. These three functions are commented and self-contained. The main role of this file is to have you ignore the supporting functions and focus on the ones that you will actually use to reproduce/modify the experiments and plots.



1) Supporting the theoretical analysis on generalized elliptical distributions

'help plot_rTME.m' for details

2) Measuring the sensitivity of TME under changes in alpha in the fixed p setting

'help alpha_sensitivity.m' for details

3) Measuring the sensitivity of TME under changes in alpha in the fixed p/n setting

'help alpha_fix_gam.m' for details

4) Investigating the performance of regularized-TME in the presence of outliers

run script 'sim_TME_outliers.m'
