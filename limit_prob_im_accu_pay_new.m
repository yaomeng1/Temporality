function prob_part = limit_prob_im_accu_pay_new(mAdj,b_array,t_array,P_coal_t,varargin)
% Computes the limit time average probability for a graph given by weighted 
% (replacement) adjacency
% matrix mAdj. 
% A first optional argument represents a distinct interaction matrix mInt
% A second optional argument can be used to pass a precomputed matrix of
% remeeting times. This is useful for analyzing different interaction
% matrices with the same replacement matrix

% parameter: b_array: column vector of b
% parameter: t_array: row vector of time t 
% return: prob: a matrix of row: prob through t; column: prob through b

c=1;
% delta=0.01;
n = length(mAdj);
w = sum(mAdj);

% pi = (w/W).';
pi = ((w+1)/ sum(w+1))';

c_term1 = zeros(1,length(t_array));
c_term2 = zeros(1,length(t_array));
c_term3 = zeros(1,length(t_array));
b_term1 = zeros(1,length(t_array));
b_term2 = zeros(1,length(t_array));
b_term3 = zeros(1,length(t_array));
b_term4 = zeros(1,length(t_array));
b_term5 = zeros(1,length(t_array));

P10 = normalizedLaplacian(mAdj)+speye(n);


Rem_t = zeros(n,n);

for step = 1:length(t_array)
 
%     Rem_t = Rem.*(1 - exp(-t_array(step)./(n/2*Rem))); % limit remeet time
%     Rem_t = sum(P_coal_t(:,:,1:step-1),3);
    Rem_t = Rem_t - diag(diag(Rem_t));
    
    c_term1(step) = (pi'.*(w./(w+1)).^2)*sum(P10.*Rem_t,2);
    c_term2(step) = (pi'.*(w./(w+1)).^2)*sum(P10.*((P10.*w)*Rem_t),2);
    c_term3(step) = (w./(w+1).^2.*pi')*sum(P10.*w.*Rem_t,2);
    
    b_term1(step) = -(pi'.*(w./(w+1)))*sum(P10.*(w'.*sum(P10.*Rem_t,2))',2);
    b_term2(step) = (pi'.*(w./(w+1)).^2)*sum(P10.*((P10.*w)*(P10*Rem_t)),2);
    b_term3(step) = (pi'.*(w./(w+1)).^2)*sum(P10.*(P10*Rem_t),2);
    
    b_term4(step) = (w./(w+1).^2.*pi')*sum(P10.*w.*(P10*Rem_t)',2);
    b_term5(step) = -((w./(w+1)).^2.*pi')*sum(P10.*Rem_t,2);
    

    Rem_t = Rem_t+ P_coal_t(:,:,step);


end


% prob_part=k/n+delta*k*(n-k)/(2*n*(n-1))*(-c*T20_t+b_array*(T21_t-T01_t));
prob_part=-c*(c_term1+c_term2+c_term3)+b_array*(b_term1+b_term2+b_term3+b_term4+b_term5);

end

