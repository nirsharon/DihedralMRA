function [groundtruth, sync_measurements] = sync_data_for_testing(n, N)
% generate N ground truth elements of D_2n and 
% a measurements of matrix consists of the blocks:
%   g_ij =  g_i g_j^{-1} .
%
% NS, Aptil 21

% ====> FULL GRAPH CASE <====

d = 2;    % embedded in O(2)

% generate random ground truth
groundtruth = rand_D2n_element(n, N);

% prepare the sync data
sync_measurements = zeros(2*N);

% The measures are from a compact group so the matrix is symmetric
for i=1:N
    ind_i = (d*i-1):(d*i);
    for j=(i+1):N
        ind_j = (d*j-1):(d*j);
        sync_measurements(ind_i,ind_j) = groundtruth(:,:,i)*(groundtruth(:,:,j)');
    end
end

% summary
sync_measurements = sync_measurements + sync_measurements.' + eye(2*N);

end

