% This script shows the workflow to fit an ACMTF model to jointly analyze two tensors (real/simulated meal challenge metabolomics data) 
%
%   L. Li, H. Hoefsloot, ..., M. A. Rasmussen, A. K. Smilde and E. Acar, Longitudinal metabolomics data analysis 
%   informed by mechanistic models, 2024 

%
% The script assumes that the data is stored as a dataset object
% (https://eigenvector.com/software/dataset-object/)
%
% Needed toolboxes:
% Tensor toolbox (https://www.tensortoolbox.org/)
% Poblano toolbox (https://github.com/sandialabs/poblano_toolbox)
% CMTF Toolbox Available on https://github.com/eacarat/CMTF_Toolbox
% Auxiliary functions: under folder ./functions

% Simulated data can be downloaded from
% https://github.com/Lu-source/project-of-challenge-test-data/tree/main/simulated_datasets/Betacell_dysfunction/Simu_6meta_8time_alpha02_betacell_balance.mat


%% load real data
load('data_real.mat')
% real data is of size 299 subjects by 8 time points by 252 features (Nightingale NMR
% measurements + insulin + cpeptide)
% Select the six metabolites (both in real and simulated data)
X=NMR;
id_sim =[251 71 73 72 61 76];
X  = X(:,:,id_sim);

% select gender: males (id_gender = 2) or females (id_gender = 1)
id_gender = find(X.class{1,1}==2);
YY        = X(id_gender,:,:); %1:females, 2:males
Meta      = Metainfo(id_gender,:); % additional information available about the subjects

% Take T0-corrected data
YYdiff   = take_diff(YY); %T0-corrected real data

% remove outliers
X = YYdiff;
outlier  = [142 79 342]; %male
%outlier = [90 335 250]; %female
[~,id,~] = intersect(str2num(X.label{1,1}), outlier);
inc      = setdiff(1:size(X,1),id);
Xreal   = X(inc,:,:);
Meta     = Meta(inc,:);

% merge BMI groups into two groups
cc = Xreal.class{1,2};
cc(find(Xreal.class{1,2}==2 | Xreal.class{1,2}==3))=2; %obese and overweight
cc(find(Xreal.class{1,2}==1 | Xreal.class{1,2}==4))=1; %underweight and normal
Xreal.class{1,11}=cc;
Xreal.classname{1,11}='BMI 2-class';

clearvars -except Xreal Meta 


%% load simulated data
load('data_simulated.mat')
% simulated data is of size 100 subjects (50 control and 50 pathology subjects) by 6 metabolites by 8 time points by 252 features 
% Select the control subjects for joint analysis
normal = find(X_orig.class{1,1}==1);
Xorig  = X_orig(normal,:,:);
Xp     = permute(Xorig, [1 3 2]); % permute the simulated to have same arrangement of modes as real data

Xsim = take_diff(Xp); %T0-corrected simulated data



%% set model parameters
R = 3;
nb_starts=50;
lambda = 0.01; % regularization coefficient (the data is normalized when fitting regularized CP models)

%% Joint analysis of real and simulated data using ACMTF
data = fit_acmtf_simreal(Xreal, Xsim, R, nb_starts);

%% Analysis of real data using the regularized CP model
[data_cp, fit_cp, out_cp] = fit_cp_ridge_real(Xreal, R, nb_starts, lambda);

%% save the results 
eval(strcat('save ACMTF_CPReal_R', num2str(R), '.mat'))


%% Example showing the weight of each component of ACMTF
server_flag=0;
legd       = {'\lambda (Real)','\sigma (Simulated)'};
f=figure;
[Fac_aligned, T1, T2] = show_spread(R, data.Fac_sorted, data.f_sorted, server_flag, legd);
set(gca,'Fontsize',15)


            
%% Example showing how to check the results, e.g., looking at correation between subject scores with meta variables 
S_ACMTF = data.Zhat{1}.U{2}; % real subjects factor of R-component ACMTF model
S_CP    = data_cp.U{1}; % real subjects factor of R-component ACMTF model
C_ACMTF = best_corr_compute(S_ACMTF,Meta);
C_CP    = best_corr_compute(S_CP,Meta);
C=[C_CP, C_ACMTF]; % correlations with meta variables
feature_label ={'HOMAIR', 'MuscleFatRatio','FatPercent', 'MuscleMass', 'Weight', 'BMI', 'Waist', 'WaistHeightRatio', 'FatMass','FatMassIndex', 'FFMI'};
figure
bar(C);
ylabel('Pearson corr. coef.')
set(gca, 'XTick',1:1:11,'XTickLabel',feature_label,'YGrid','on')
hold on; xtickangle(45)
lng={'CP', 'ACMTF'};
legend(lng)
set(gca,'Fontsize',15)


%% Replicability for different number of components
for r = 1:4
    [FMS_acmtf_real{r},FMS_acmtf_sim{r},FMS_cp_real{r}]  = script_simreal_replicability(Xreal, Xsim, r, nb_starts, lambda);
end
plot_fms_replicability(FMS_cp_real)

