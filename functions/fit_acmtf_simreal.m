
% This function is used to fit the real and simulated data using an ACMTF model

function data = fit_acmtf_simreal(Xreal, Xsim, R, nb_starts)

%centering and scaling
Xrealp = preprocess_centerscale(Xreal,1,1);
Xsimp  = preprocess_centerscale(Xsim, 1,1 );

%form coupled object
X{1} = permute(Xrealp, [3 1 2]);
X{2} = permute(Xsimp, [3 1 2]);

X1 = tensor(X{1}.data);
W1 = tensor(~isnan(X{1}.data));
X1(find(W1==0))=0;
Z.object{1} = X1;

X2 = tensor(X{2}.data);
W2 = tensor(~isnan(X{2}.data));
Z.object{2} = X2;

Z.modes = {[1 2 3]; [1 4 5]};
Z.size = [size(X{1}) size(X{2},2) size(X{2},3)];
Z.miss{1} = W1;
Z.miss{2} = W2;
for i=1:2
    Z.object{i}=Z.object{i}/norm(Z.object{i});
end

%% fit ACMTF
options = ncg('defaults');
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-8;
options.RelFuncTol   = 1e-8;
 
beta   = [0.001 0.001]; % l1-penalty parameter for the higher-order tensors

for i = 1:nb_starts
  [Fac{i},~,out{i}] = acmtf_opt(Z,R,'init','random','alg_options',options, 'beta', beta, 'alg', 'ncg');    
end
%%
P =2;
for i=1:nb_starts
    f(i) = out{i}.F;
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

data.Fac_sorted=Fac_sorted; 
data.f_sorted=ff; 

Zhat = Fac_sorted{1};
l_rec = zeros(P, R);
for p = 1:P
    temp        = normalize(Zhat{p});
    l_rec(p,:)  = temp.lambda;    
end
for p=1:P
    data.Zhat{p}= normalize(Zhat{p});
    data.fit(p) = 100- (norm(Z.miss{p}.*tensor(full(data.Zhat{p})-Z.object{p}))^2/norm(Z.miss{p}.*Z.object{p})^2*100);
end
for p=1:P
    for r=1:R
        [iii,jjj] = max(abs(data.Zhat{p}.U{3}(:,r)));
        if data.Zhat{p}.U{3}(jjj,r)<0
            data.Zhat{p}.U{3}(:,r)=- data.Zhat{p}.U{3}(:,r);
            data.Zhat{p}.U{1}(:,r)=- data.Zhat{p}.U{1}(:,r);            
        end
    end
end

