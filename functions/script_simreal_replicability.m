function [FMS_acmtf_real, FMS_acmtf_sim, FMS_cp_real] = script_simreal_replicability(Xreal, Xsim, R, nb_starts, lambda)

% This function checks the replicability of R-component ACMTF/CP model by the following steps
%   Step1: form a random ten-fold of the samples in the first mode of data tensor Xreal 
%   Step2: leave one fold of the data out 
%   Step3: fit an R-component ACMTF/CP model to the remaining 90% - repeat this by leaving out each fold once
%   Step4: compute the factor match scores between the best runs of the ACMTF/CP models (in the second and third modes for the real data part
%          , but all modes in the simulated data part). 
%   Step5: Repeat Step1-4 ten times
% Return in total 450 FMS scores.

%% Reorder subjects
sub_low    = find(Xreal.class{1,2}==1 | Xreal.class{1,2}==4); % Lower BMI subjects
sub_high   = find(Xreal.class{1,2}==2 | Xreal.class{1,2}==3); % Higher BMI subjects
index_perm = [sub_low,sub_high];
Xreal_perm          = Xreal(index_perm,:,:);


%% form splits and fit ACMTF/CP
nb_splits      = 10;
repeats        = 10;
FMS_acmtf_real = [];
FMS_acmtf_sim  = [];
FMS_cp_real    = [];
for r = 1:repeats
    index_perm = [randperm(length(sub_low)),randperm(length(sub_high))+length(sub_low)];
    for s = 1:nb_splits
        index_rem = index_perm(s:nb_splits:end); %randomly remove 1/nb_splits subjects (equally many in each group)
        % data with remaining subjects
        Xreal_inc{r}{s}  = Xreal_perm(setdiff(index_perm, index_rem),:,:);     
        data_ACMTF{r}{s} = fit_acmtf_simreal(Xreal_inc{r}{s}, Xsim, R, nb_starts);
        index_uni(s) = check_uni_coupledtensor(data_ACMTF{r}{s},R);
        [Fac_CP{r}{s}, fit_cp{r}{s}] = fit_cp_ridge_real(Xreal_inc{r}{s}, R, nb_starts,lambda);
    end
    % select unique Facs in ACMTF
    index_sel=find(index_uni==1);
    for j=1:length(index_sel)
        Fac_ACMTF_real{r}{j}=permute(data_ACMTF{r}{index_sel(j)}.Zhat{1},[2,3,1]);
        Fac_ACMTF_sim{r}{j}=permute(data_ACMTF{r}{index_sel(j)}.Zhat{2},[2,3,1]);
    end
    
    % compute FMS
    fms_1{r}    = compute_pairwise(Fac_ACMTF_real{r});
    fms_2{r}    = compute_pairwise_v1(Fac_ACMTF_sim{r});
    fms_cp{r}   = compute_pairwise(Fac_CP{r});
    FMS_acmtf_real = [FMS_acmtf_real;fms_1{r}];
    FMS_acmtf_sim  = [FMS_acmtf_sim;fms_2{r}];
    FMS_cp_real    = [FMS_cp_real;fms_cp{r}];
end



eval(strcat('save FMS_replicability_ACMTF_CP_R', num2str(R), '.mat'))



%%
% This function is used to compute FMS between two pairs of factors using second and third modes, i.e.,
% metabolites and time modes, to check the replicability of factors extracted from real data by ACMTF model.

function fms = compute_pairwise(Fac) 

nb_splits = length(Fac);
for i=1:nb_splits
    for j=i+1:nb_splits
        Sim(i,j)= score(ktensor(Fac{i}.U([2:3])), ktensor(Fac{j}.U([2:3])));
    end
end
fms = Sim(find(Sim>0));



% This function is used to compute FMS between two pairs of factors using all modes, i.e.,
% subjects, metabolites and time modes, to check the replicability of factors extracted 
% from simulated data by ACMTF model.

function fms = compute_pairwise_v1(Fac)

nb_runs = length(Fac);
for i=1:nb_runs
    for j=i+1:nb_runs
        Sim(i,j)= score(ktensor(Fac{i}.U), ktensor(Fac{j}.U));
    end
end
fms = Sim(find(Sim>0));
