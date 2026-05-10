# Code and Data for OCRL

Overview

1. Directory `Simulation` contains the main functions required to perform simulations for the proposed online change rate learning (OCRL), as well as several classical offline and online comparative methods. Specifically,

1.1 The sub-directory `OCRL` contains the functions for the OCRL algorithm.

       The R script `basic functions.R` provides the essential functions for both OCRL and the classical local polynomial regression.

       The R scripts `Mean1-st_dense_OCRL.R` and `Mean1-st_sparse_OCRL.R` estimate the first-order derivative of the mean function under dense and sparse settings.

       The R scripts `Cov-1st_C_dense_OCRL.R`, `Cov-1st_C_sparse_OCRL.R` and `Cov1-st_dense_sparse_OCRL.R` estimate the first-order derivative of the covariance function under the dense and sparse settings.

1.2 The sub-directory `Kernel` contains the functions for the classical local polynomial regression. 

      The R scripts `Mean1-st_dense_kernel.R`, `Mean1-st_sparse_kernel.R`, `Cov-1st_C_dense_kernel`, `Cov-1st_C_sparse_kernel.R` and `Cov1-st_dense_sparse_kernel.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

1.3 The sub-directory `SGD` contains the functions for the stochastic gradient descent method.

      The R scripts  `SGDo_Mean1-st.R` and `SGDo_Cov1-st.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

1.4 The sub-directory `Spline` contains the functions for conducting simulations using the spline-based method.

      The R scripts `Spline_Mean1-st_fun.R`  and `Spline_Cov1-st_fun.R` provide the essential functions for the spline-based method.

      The R scripts `Spline_Mean1-st.R` and `Spline_Cov1-st_fun.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

1.5 The sub-directory `Wavelet` contains the functions for conducting simulations using the wavelet-based method.
      
      The R script `Waveopt_Cov1-st.R` provides the essential functions for the wavelet-based method.

      The R scripts `Waveopt_Mean1-st.R` and `Waveopt_Cov1-st_fun.R` estimate the first-order derivatives of the mean and covariance functions under the dense and sparse settings.

2. Directory `Realdata` contains the preprocess for the two real data examples which are available in the following links:

- ***Power_Consumption*** : https://www.eia.gov/electricity.
- ***US_Traffic_Accident*** :  https://smoosavi.org/datasets/us_accidents.

      The R scripts `dense_mean1-st_online_and_kernel.R`and `sparse_mean1-st_online_and_kernel.R` provide the codes for analyzing the hourly power consumption data in megawatts.
      
      The R scripts `count_USAcc_dense_mean1-st_online_and_kernel.R` and `count_USAcc_sparse_mean1-st_online_and_kernel.R` provide the codes for analyzing the traffic accident data from the United States.

# Reproducibility Workflow

The code was developed under R version 4.4.2 and depends only on widely used CRAN packages, including doParallel, dplyr, fda, fda.usc, foreach, KernSmooth, Matrix, mgcv, pracma, and wavethresh.

To reproduce the numerical results in the manuscript, users may follow the workflow below.

1. Simulation Studies
(1). Navigate to the directory Simulation.
(2). Choose one of the above subdirectories depending on the method of interest.
(3). Run the corresponding scripts according to the sampling scheme:
*_dense_*.R for dense functional data, *_sparse_*.R for sparse functional data,*_Mean1-st_*.R for the first-order derivative estimation of the mean function, *_Cov1-st_*.R for the first-order partial derivative of the covariance function.

The following example reproduces the dense simulation result for the first-order derivative of the mean function using the proposed OCRL method.

setwd("Simulation/OCRL")
source("basic functions.R")
source("Mean1-st_dense_OCRL.R")

Running the above scripts generates the estimated results and the corresponding performance measures reported in the simulation studies.

2. Real Data Applications

The real data analyses are located in the directory Realdata. The empirical analyses can be reproduced by running the R.file in the sub-directory.





