
% This function is used to fit the real  data using an regularized CP model
% lambda is the regularization coefficient
% we consider Tikhonov regularization

function [data, fit, out_flag] = fit_cp_ridge_real(Xreal, R, nb_starts,lambda)

X = preprocess_centerscale(Xreal,1,1);

X = tensor(X.data);
W = tensor(~isnan(X.data));
X(find(W==0)) = 0;

X=X/norm(X);


%% fit CP with ridge
options = ncg('defaults');
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-8;
options.RelFuncTol   = 1e-8;
 
for i = 1:nb_starts
  [Fac{i},~,out{i}] = cp_wopt_reg(X,W, R,lambda,'init','randn','opt_options',options, 'opt', 'ncg');    
end
%%
for i=1:nb_starts
    f(i) = out{i}.f;
end
[ff, index] = sort(f,'ascend');
%check the stopping condition
if (out{index(1)}.ExitFlag==0) | (out{index(1)}.ExitFlag==3)
   flag_stop = true;
else 
   flag_stop = false;
end
Fac_sorted = Fac(index);
out_sorted = out(index);
out_flag.ExitFlag=out{index(1)}.ExitFlag;
out_flag.fms=score(Fac_sorted{1},Fac_sorted{2},'lambda_penalty',false);
out_flag.ff=ff;
out_flag.Fac_sorted=Fac_sorted;
out_flag.out_sorted=out_sorted;
data = normalize(Fac_sorted{1});
fit = 100- (norm(W.*tensor(full(data)-X))^2/norm(W.*X)^2*100);
