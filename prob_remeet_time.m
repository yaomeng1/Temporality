function prob_time_mat = prob_remeet_time(mAdj,time_length)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% IM update

% time_length = 1500;
n = length(mAdj);
mAdj = mAdj+ eye(n);

P10 = normalizedLaplacian(mAdj)+speye(n);

% 
% prob_time_mat = zeros(n,n,time_length);
% 
% prob_time_mat(:,:,1) = 0.5*P10+0.5*P10';
% 
% 
% for i = 2:time_length
%     prob_slice_t = prob_time_mat(:,:,i-1);
%     for j = 1:n
%         for k = 1:n
%             minus = zeros(1,n);
%             minus(1,k) = P10(j,k);
%             prob_time_mat(j,k,i)=(P10(j,:)-minus)*prob_slice_t(:,k);
%         end
%     end
%      prob_time_mat(:,:,i) = 0.5*(prob_time_mat(:,:,i)'+prob_time_mat(:,:,i));
% end 

Prob_mat = zeros(n,n,time_length);
Prob_mat_0 =  eye(n); %% use 0 step to meet
Prob_mat(:,:,1) = 1/2*P10*Prob_mat_0+1/2*Prob_mat_0*P10';
Prob_mat(:,:,1) = Prob_mat(:,:,1) - diag(diag(Prob_mat(:,:,1)));
for i=2:time_length
    Prob_mat(:,:,i) = 1/2*P10*Prob_mat(:,:,i-1)+1/2*Prob_mat(:,:,i-1)*P10';
    Prob_mat(:,:,i) = Prob_mat(:,:,i) - diag(diag(Prob_mat(:,:,i)));
end

prob_time_mat = Prob_mat;
end

