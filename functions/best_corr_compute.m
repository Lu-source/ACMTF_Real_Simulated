% This auxiliary function is used to compute the correlations achieved by 
% using the component of a model, which gives the best correlations.

function [C,index_best_r] = best_corr_compute(Sub_Fac,Meta)


meta_feature = [11:14 25:31];
Meta = Meta(:,meta_feature);
R = size(Sub_Fac,2);

for i = 1:size(Meta,2)
    idnan = isnan(Meta.data(:,i));
    for k = 1:R
        if sum(idnan)<50
            [C(i,k), ~] = corr(Sub_Fac(find(idnan==0),k),Meta.data(find(idnan==0),i));
        else
            C(i,k)  = NaN;
        end
    end
end

[~,index_best_r_temp] = max(abs(C),[],2);
index_best_r = mode(index_best_r_temp);
   
    if C(6,index_best_r)<0
        
          C=-C;
    end
    
    C=C(:,index_best_r); % best correlation

end



