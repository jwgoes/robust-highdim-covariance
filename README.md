# Code for "Robust Sparse Covariance Estimation by Thresholding Tyler’s M-Estimator"
Code for "Robust Sparse Covariance Estimation by Thresholding Tyler’s M-Estimator", by John Goes, Gilad Lerman and Boaz Nadler (https://arxiv.org/pdf/1706.08020.pdf).

*******************************************************************************************************************

Paper abstract: Estimating a high-dimensional sparse covariance matrix from a limited number of samples is a fundamental problem in contemporary data analysis. Most proposals to date, however, are not robust to outliers or heavy tails. Towards bridging this gap, in this work we consider estimating a sparse shape matrix from n samples following a possibly heavy tailed elliptical distribution. We propose estimators based on thresholding either Tyler’s M-estimator or its regularized variant. We derive bounds on the difference in spectral norm between our estimators and the shape matrix in the joint limit as the dimension p and sample size n tend to infinity with p/n → γ > 0. These bounds are minimax rate-optimal. Results on simulated data support our theoretical analysis.

*******************************************************************************************************************


This page describes four matlab functions that can be used to replicate experiments in the above paper:

1) Supporting the theoretical analysis (produces panel of Figure 1 in paper):

'help plot_rTME_exp.m' for details

2) Measuring the sensitivity of TME under changes in alpha in the fixed p/n setting (Figure 2 in paper):

'help many_alphas_exp.m' for details

3) Measuring the sensitivity of TME under changes in alpha in the fixed p setting (Figure 3 in paper):

'help alpha_sensitivity_exp.m' for details

4) Investigating the performance of regularized-TME in the presence of outliers (Figure 4 in paper):

run script 'sim_TME_outliers_exp.m'

