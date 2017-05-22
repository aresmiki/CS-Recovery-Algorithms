========================================================================
       			CS Recovery Algorithms
========================================================================
CS Recovery Algorithms Toolbox v0.1
Authors: He Liu & Lin Yan (https://github.com/aresmiki/CS-Recovery-Algorithms.git)

Introduction
--------------

This toolbox implements several CS recovery Algorithm
methods, as described in [1,2,3,4,5,6,7,8,9,10,11,12,13,14].  These include:

1. Orthogonal Matching Pursuit
2. Compressive Sampling Matching Pursuit
3. Fast Iterative Shrinkage-Thresholding Algorithm
4. Iterative Hard Thresholding algorithms for compressive sensing
5. Iteratively Reweighted Least Square
6. Iterative Shrinkage-Thresholding Algorithm
7. Null-Space Reweigthted Approximate l0-Pseudonorm Algorithm
8. Reweighted L1 Minimization Algorithm
9. Robust Smoothed l0-Pseudonorm Algorithm
10. L1_SplitBregmanIteration
11. Smoothed l0-Pseudonorm Algorithm
12. Minimization of Approximate Lp Pseudonorm Using a Quasi-Newton Algorithm

Installation
-------------

0. You will need the MATLAB L1-MAGIC toolbox for solving the convex 
   optimization programs central to compressive sampling.

1. Put the contents of the CS Recovery Algorithms toolbox somewhere (say,
   $HOME/matlab/).

2. Add the new directories to your path permanently; e.g., add the
   following to your startup.m:
     addpath ~/matlab/CS Recovery Algorithms;

Basic command reference
------------------------
1. CS_OMP
2. CS_CoSaMP
3. CS_FISTA
4. CS_IHT
5. CS_IRLS
6. CS_ISTA
7. CS_NSRAL0
8. CS_RL1
9. CS_RSL0
10. CS_SBIL1
11. CS_SL0
12. CS_UALP


References
-------------
[1] Joel A. Tropp and Anna C. Gilbert Signal Recovery From Random Measurements Via Orthogonal Matching Pursuit,IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12.
[2] Needell D,Tropp J A CoSaMP:Iterative signal recovery from incomplete and inaccurate samples[J].Applied and Computation Harmonic Analysis,2009,26:301-321.
[3] D.Needell, J.A. Tropp.CoSaMP:Iterative signal recoveryfrom incomplete and inaccurate samples[J].Communications of theACM,2010,53(12):93-100.
[4] A. Beck and M. Teboulle, "A fast iterative shrinkage-thresholding algorithm for linear inverse problems," SIAM J. Imaging Sciences, vol. 2, no. 1, pp. 183-202, 2009.
[5] Blumensath T, Davies M E. Iterative hard thresholding for compressed sensing[J]. Applied & Computational Harmonic Analysis, 2009, 27(3):265-274.
[6] Blumensath T, Davies M E. Iterative Thresholding for Sparse Approximations[J]. Journal of Fourier Analysis and Applications, 2008, 14(5):629-654.
[7] Chartrand and W. Yin, "Iteratively Reweighted Algorithms for Compressed Sensing," 2008.
[8] I. Daubechies, M. Defrise, and C. D. Mol, "An iterative thresholding algorithm for linear inverse problems with a sparsity constraint," Comm. Pure Appl. Math., vol. 57, pp. 1413-1457, 2004.
[9] J. K. Pant, W.-S. Lu, and A. Antoniou,"Reconstruction of sparse signals by minimizing a re-weighted approximate l0-norm in the null space of the measurement matrix," IEEE Inter. Midwest Symp. on Circuits-Syst, pp. 430-433, 2010.
[10] Cand¨¨s E J, Wakin M B, Boyd S P. Enhancing sparsity by reweighted L1 minimization.[J]. Journal of Fourier Analysis & Applications, 2007, 14(5):877-905.
[11] H. Mohimani, M. Babie-Zadeh, and C. Jutten,"A fast approach for overcomplete sparse decomposition based on smoothed l0-norm," IEEE Trans. Signal Process., vol. 57, no. 1, pp. 289-301, Jan. 2009.
[12] Yin W, Osher S, Goldfarb D, et al.Bregman Iterative Algorithms for L1 Minimization with Applications to Compressed Sensing[J]. Siam Journal on Imaging Sciences, 2008, 1(1):143-168.
[13] H. Mohimani, M. Babie-Zadeh, and C. Jutten,"A fast approach for overcomplete sparse decomposition based on smoothed l0-norm," IEEE Trans. Signal Process., vol. 57, no. 1, pp. 289-301, Jan. 2009.
[14] Pant J K, Lu W S, Antoniou A.Unconstrained regularized Lp -norm based algorithm for the reconstruction of sparse signals[C] IEEE International Symposium on Circuits and Systems. IEEE, 2011:1740-1743.
