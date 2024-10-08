function  [Fac_aligned, T1, T2] = show_spread(R, Fac_f100, ff, server_flag, legd)

min_val = min(ff);
nb_rep  = length(find((ff-min_val)<=1e-6));
%nb_rep =5;
P = length(Fac_f100{1});
for i=1:nb_rep
    c   = abs(corr(Fac_f100{1}{1}.U{1},Fac_f100{i}{1}.U{1}));
    ind = zeros(R,1);
    for k=1:R
        temp = find(c(k,:)>0.95);
        if isempty(temp)
            ind(k) = 0;
        else
            ind(k) = temp(1);
        end           
    end
    rm_clm = setdiff(1:R,setdiff(ind, 0));
    if length(rm_clm)>1
        ind = place_remaining(c,ind,R);
    else
        ind(find(ind==0))=setdiff(1:R,setdiff(ind, 0)); 
    end
    Fac_aligned{i} = perm_match(Fac_f100{i},ind, P);    
    clear ind;
end


for i=1:nb_rep
    t1 = normalize(Fac_aligned{i}{1});
    t2 = normalize(Fac_aligned{i}{2});
    T1(i,:) = t1.lambda';
    T2(i,:) = t2.lambda';   
    if P==3
        t3 = normalize(Fac_aligned{i}{3});
        T3(i,:) = t3.lambda';   
    end
end

if ~server_flag
    if nb_rep>1
        if P==2
            x1 = 0.85:1:R-0.15;
            x2 = 1.15:1:R+0.15;
            lims = [0 R+1 0 1];
            %[rmin, rmax] = bar_wrange(T1,T2,{'\lambda(ERP)','\sigma(fMRI)'},x1,x2,lims);
            [rmin, rmax] = bar_wrange(T1,T2,legd,x1,x2,lims);
           
            %set(gca,'Ygrid','on', 'YTick', 0.1:0.2:0.7);
            set(gca,'Ygrid','on', 'YTick', 0.1:0.2:1);
        elseif P==3
            x1 = 0.775:1:(R-0.225); 
            x2 = 1:1:R;
            x3 = 1.225:1:R+0.225;
            lims = [0.5 R+1 0 0.95];
            [rmin, rmax] = bar_wrange_multiple(T1,T2,T3,legd, x1, x2, x3,lims);
    
            %[rmin, rmax] = bar_wrange_multiple(T1,T2,T3,{'\lambda (ERP)','\sigma (fMRI)','\gamma (sMRI)'}, x1, x2, x3,lims);
            set(gca,'Ygrid','on',  'YTick', 0.1:0.2:1);
        end
        xlabel('Components');    
    end
end
if P==3
    varargout=T3;
end


function Y = perm_match(X, ind,P)

for k=1:length(X{1}.U)
    Y{1}.U{k} = X{1}.U{k}(:,ind);
end
Y{1}.lambda = X{1}.lambda(ind);
Y{1} = ktensor(Y{1}.lambda, Y{1}.U);

for k=1:length(X{2}.U)
    Y{2}.U{k} = X{2}.U{k}(:,ind);
end
Y{2}.lambda = X{2}.lambda(ind);
Y{2} = ktensor(Y{2}.lambda, Y{2}.U);

if P==3
    for k=1:length(X{3}.U)
        Y{3}.U{k} = X{3}.U{k}(:,ind);
    end
    Y{3}.lambda = X{3}.lambda(ind);
    Y{3} = ktensor(Y{3}.lambda, Y{3}.U);
end

function  ind = place_remaining(c,ind,R)

for i=1:R
    if ind(i)==0
        rm_clm = setdiff(1:R,setdiff(ind, 0));
       [ii,jj]= max(c(i,rm_clm));
       ind(i) = rm_clm(jj);
       clear ii jj
    end
end
