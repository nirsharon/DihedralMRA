function[ g_array ] = sync_D2n(H, n)
% 
% The eigenvectors synchronization algorithm -- adapted to D_2n
%
% N.S, April 21

% the affinity matrix
d = 2;                  % 2 X 2 matrix representation for D_2n
N = size(H,1)/d;        % number of unknown elements

% if H is not symmetric, choose the upper part of H
if norm(H-H')>1e-14
    for j = 1:N
        ind_j = (d*j-1):(d*j);
         H(ind_j,ind_j) = eye(d);
        for l=(j+1):N
                ind_l = (d*l-1):(d*l);
                H(ind_l,ind_j) = H(ind_j,ind_l)';
        end
    end
end


%extract eigenvectors
try
    [vecs, ~] = eigs(H,d);
catch
    [vecs, ~] = eigs(H,d); %%% to avoid troubles
end

% rounding -- first to rotation
rotations_arr = zeros(d,d,N);
for j=1:N
    ind1 = 1+(j-1)*d;        
    B = vecs(ind1:(ind1+d-1),:);
    [u , ~, v] = svd(B);
    rotations_arr(:,:,j) = D2n_rounding(u*v', n);
    %-----
%     ind1 = 1+(j-1)*d;        
%     vv_so(ind1:(ind1+d-1),:) = u*v';
%     rot_arr(:,:,j) = u*v';
    %-----
    B2 = vecs(ind1:(ind1+d-1),[2,1]);
    [u2 , ~, v2] = svd(B2);
    rotations_arr2(:,:,j) = D2n_rounding(u2*v2', n);
end
%g_array = rotations_arr;

% 
% % debugging
% %Q = rot_arr(:,:,1);
for j=1:N
% %    rot_arr(:,:,j) = rot_arr(:,:,j)*Q';
% %    gg_arr(:,:,j) =  D2n_rounding(rot_arr(:,:,j), n);
    ind1 = 1+(j-1)*d;        
    vv(ind1:(ind1+d-1),:) = rotations_arr(:,:,j);
    vv2(ind1:(ind1+d-1),:) = rotations_arr2(:,:,j);
 end
remainder1 = norm(H*vv-N*vv);
remainder2 = norm(H*vv2-N*vv2);


%norm(H*vv_so-N*vv_so)
% rounding -- to D_2n version 1 (reflect)
% for j=1:N
%     rotations_arr(:,:,j) = [rotations_arr(:,2,j), rotations_arr(:,1,j)];
% end
% g_array2 = round_rotation_to_D2n(rotations_arr2, n);
% for j=1:N
%    ind1 = 1+(j-1)*d;        
%    vv3(ind1:(ind1+d-1),:) = g_array1(:,:,j);
% end
% norm(H*vv3-10*vv3)
% 
% remainder1 = grouped_proc_cost_grad(g_array1, kron(ones(N),eye(2)), H, ones(2*N));
% remainder2 = grouped_proc_cost_grad(g_array2, kron(ones(N),eye(2)), H, ones(2*N));
% 
if remainder1<remainder2
    g_array = rotations_arr;
else
    g_array = rotations_arr2;
end

%g_array  = gg_arr;
end

