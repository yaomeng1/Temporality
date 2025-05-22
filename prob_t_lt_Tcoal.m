function P_coal_t = prob_t_lt_Tcoal(prob_time_mat)
% the intergral of probability density function from t to infinity,
% that is also P[t<Tcoal]
% prob_time_mat: (i,j,t) prob density for any i and j at time t (from t=1)

prob_size = size(prob_time_mat);
n = prob_size(1);
t_length = prob_size(3);
% prob_mat = zeros(n,n,1+t_length);
% prob_mat(:,:,2:end) = prob_time_mat;
% 
% P_coal_t = zeros(n,n,1+t_length); % integral from t to infinity
% 
% 
% P_coal_t(:,:,end) = prob_mat(:,:,end);
% for t = 1:t_length
%     P_coal_t(:,:,end-t) = P_coal_t(:,:,end-t+1)+prob_mat(:,:,end-t);
% end

P_coal_t = zeros(n,n,1+t_length); % integral from t to infinity
P_coal_t(:,:,end) = 0;
for t = 1:t_length
    P_coal_t(:,:,end-t) = P_coal_t(:,:,end-t+1)+prob_time_mat(:,:,end-t+1);
end

end
