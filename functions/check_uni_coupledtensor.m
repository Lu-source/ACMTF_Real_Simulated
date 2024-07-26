function [index, index_1] = check_uni_coupledtensor(data,R)
% This function is used to numerically check whether an ACMTF model 
% produces unique factorizations or not

Fac_aligned = show_spread(R, data.Fac_sorted, data.f_sorted,1);
if length(Fac_aligned)<2
    disp ('Need more starts')
    index=0;
else
    ll=0;
    for jj = 1:length(Fac_aligned)
        for kk = jj+1:length(Fac_aligned)
            ll=ll+1;
            FMS_score_tensor1(ll) = score(Fac_aligned{jj}{1},Fac_aligned{kk}{1},'lambda_penalty',false);
            FMS_score_tensor2(ll) = score(Fac_aligned{jj}{2},Fac_aligned{kk}{2},'lambda_penalty',false);
        end
    end
    if min(FMS_score_tensor1)>=0.95
        index_1 = 1; % factorization with unique factors--tensor1
    end
    if min(FMS_score_tensor2)>=0.95
        index_2 = 1; % factorization with unique factors--tensor2
    end
    if  index_1 == 1 &&  index_2 ==1
        index=1;
    else
        index=0;
    end
end