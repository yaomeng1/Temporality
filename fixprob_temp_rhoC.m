% prob vector evolution
%%
clear
clc
tic 
warning('off');
global k_ini n ratio

load('re_200_k6_4_snapmatrix.mat');
limit_idx = 100;

delta=0.01;
N=200;
delta_t=1;
time_length=3000;
time_array=0:delta_t:time_length;     % time_array must bt row vector
% theory prob density and P[t<T_coal]
% [a_mat,tau_a,b_mat,tau_b] = prob_fit_2expon(prob_time_mat,P_coal_t);
% Center_tau = tau_array(P_coal_t,round(n));
%%

b_array = 7.407;     % b must be scalar 
% %b*=7.407 for g=1e3 
% %b*= 6.25 for g=2e3 
% %b*= 4.95 for g=1e4 

p_final_fix_sf = zeros(1,length(b_array));
for b_index=1:length(b_array)
    
b = b_array(b_index); % b must be scalar

% -------------begin iteration computation of a single b------------- 
iter_num=500;
prob_total = zeros(iter_num+1,N+1);
prob_total(1,2) =1; %beginning p1=1;else,pi=0
for idx=1:iter_num
   
    prob_idx =  prob_total(idx,:);  
    temp_p0_pn=[prob_idx(1),prob_idx(end)]

    %%
    if mod(idx,limit_idx)==0
        mAdj = maxComp_snapshot(matrix_snap,limit_idx);
    else
    mAdj = maxComp_snapshot(matrix_snap,mod(idx,limit_idx));
    end
    prob_time_mat = prob_remeet_time(mAdj,time_length);
    P_coal_t = prob_t_lt_Tcoal(prob_time_mat);
    
    %%
    n = length(mAdj);
    prob_vec_sub = zeros(n-1,n+1); 
    k_array=[1:n-1]';

    g = round(2*1e3/n);
    
    prob_part = limit_prob_im_accu_pay_new(mAdj,b,time_array,P_coal_t);
    mFreq=k_array/n + delta*k_array.*(n-k_array)/(2*n*(n-1))*prob_part; 
    mFreq(mFreq>1)=1;
    Center_tau = tau_array(P_coal_t,mFreq,round(n));
    fix_prob = mFreq(:,end);
    mFreq_g = mFreq(:,g);
    
    for k = 1:n-1
        k_ini=k;
%         if k<n/2
        prob_vec_sub(k,end) = fix_prob(k)-fix_prob(k)*exp(-g/Center_tau(round(n-k)));
        prob_vec_sub(k,1) = 1-fix_prob(k)-(1-fix_prob(k))*exp(-g/Center_tau(round(k)));
%         else
%         prob_vec_sub(k,end) = fix_prob(k)-fix_prob(k)*exp(-g/Center_tau(n-k));
%         prob_vec_sub(k,1) = 1-fix_prob(k)-(1-fix_prob(k))*exp(-g/Center_tau(end));
%         end 
        pn = prob_vec_sub(k,end);
        p0 = prob_vec_sub(k,1);

        if 1-p0-pn==0
            d=0;
            a=0;
        else
            ratio = n*(mFreq_g(k)-pn)/(1-p0-pn);
%             d = fzero(@root_b,0.5);
            [d,fval] = fsolve(@root_b,0.95); % initial piont is 0.95
            if d>0
            a=(1-p0-pn)*(1 - d)/(1 + d - d^k - d^(-k + n));
            for i = 2:n
            prob_vec_sub(k,i) = a*d^(abs(i-1-k));
             end
            else
                prob_vec_sub(k,2:n)=(1-p0-pn)/(n-1);
            end
        end
        


    end
%     prob_vec_sub(prob_vec_sub<0)=0;
    for i=0:N
        % the total net contain i cooperator
        pi = prob_idx(i+1); 
    if pi~=0
    for k=0:n
        if k<=i && k>=n+i-N
            % k of the cooperators are selected into subnet max component
%             p_sik = nchoosek(i,k)*nchoosek(N-i,n-k)/nchoosek(N,n);
            p_sik = p_ikApprox(N,n,i,k);
            % vector(from k coop start, to 0-n prob): probVec_k_g
            if k==0 
                probVec_k_g=[1,zeros(1,n)];
            elseif k==n
                probVec_k_g=[zeros(1,n),1];
            else
                probVec_k_g = prob_vec_sub(k,:);
            end
            prob_total(idx+1,i-k+1:i-k+1+n) = prob_total(idx+1,i-k+1:i-k+1+n)+pi*p_sik*probVec_k_g;   

        end
    end
    end
    end 
    if 1-prob_total(idx+1,1)-prob_total(idx+1,end)<1e-5                                                                             
        break
    end
end

prob_total = prob_total(1:idx+1,:); % get rid of zeros of unused index
% save('./prob_total_vec_re.mat','prob_total');
p_final_fix_sf(b_index) = prob_total(end,end);
end

sf_fixprob_b = [b_array;p_final_fix_sf];

% save('./re_fixprob_b_array_dt_g1e3.mat','sf_fixprob_b');
