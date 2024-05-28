The repository contains scripts showing how to jointly analyze real and simulated meal challenge metabolomics data using an ACMTF model and how to fit a CP model to real metabolomics data. The data is assumed to be arranged as a third-order tensor with modes: subjects, time and metabolites. 

While the data we have used is not available publicly, the same pipeline can be used to analyze similar time-resolved data sets. For more details about the data and the analysis, please see the full paper:
 L. Li, H. Hoefsloot, ..., M. A. Rasmussen, A. K. Smilde and E. Acar, Longitudinal metabolomics data analysis informed by mechanistic models, 2024 



The implementation makes use of the following toolboxes/packages: 
- Brett W. Bader, Tamara G. Kolda and others. MATLAB Tensor Toolbox, Version 3.1. Available on https://www.tensortoolbox.org, 2020
- Evrim Acar et al, “Structure-revealing data fusion”, BMC Bioinformatics, 15:239, 2014. CMTF Toolbox Available on https://github.com/eacarat/CMTF_Toolbox
- Daniel M. Dunlavy, Tamara G. Kolda, and Evrim Acar, “Poblano v1.0: A Matlab Toolbox for Gradient-Based Optimization”, 2010. Available online at https://github.com/sandialabs/poblano_toolbox
- Eigenvector Research, DataSet Object, available online at https://eigenvector.com/software/dataset-object/
- Auxiliary functions are under the folder 'functions'


### 'script_ACMTF_CP_real_sim.m' is the main function. 
It fits the ACMTF and CP models to the data and shows how to relate the subject factor matrix to various meta variables. It shows how to plot the weightes of each component of ACMTF model. It also calls the script_simreal_replicability.m function to assess the replicability of ACMTF and CP models.

### Folder 'functions' contains auxiliary functions
- functions/script_simreal_replicability.m checks the replicability of the factors extracted by an R-component ACMTF/CP model. 
- functions/plot_fms_replicability.m is used to plot the factor match score values from the replicability check
- functions/fit_acmtf_simreal.m is the code to fit an ACMTF model to real and simulated metabolimics data, including preprocessing steps
- functions/fit_cp_ridge_real.m is the code to fit a CP model (with Tikhonov regularization) to real simulated metabolimics data, including preprocessing steps
- functions/preprocess_centerscale.m is the code to preprocess the data
- functions/show_spread.m and functions/bar_wrange.m are auxilary function used to plot the weight of each component of ACMTF model
- subfolder 'functions/auxiliary_regularization' includes the auxiliary functions for running CP model with Tikhonov regularization




  
