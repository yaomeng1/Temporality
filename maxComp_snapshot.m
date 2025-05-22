function subnet_matrix = maxComp_snapshot(Adjmatrix,idx)

m=Adjmatrix(idx,:,:);
N = length(m);
idx_array=1:N;
idx_mat=idx_array(ones(N,1),:);
m=reshape(m,N,N);
max_comp=max(idx_mat(m>0));  % max component node
subnet_matrix = m(1:max_comp,1:max_comp);

end

