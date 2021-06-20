function [est] = MRA_D2n_sync(Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Solving MRA over dihedral group using cc+sync
%
% INPUT
% Y -- array data vectors as its columns
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[L,N] = size(Y);

% converting the group elements to mat form
base_ang = 2*pi/L;
mat = @(x) [cos(x),-sin(x); sin(x), cos(x)];

% estimating the relative group action
relative_actions = zeros(2*N);
for i=1:N
    indi = (2*(i-1)+1):(2*i);
    refi = Y(:,i);
    for j=(i+1):N
        indj = (2*(j-1)+1):(2*j);
        
        refj = Y(:,j);
        c_ij = circcorr(refi, refj);
        [c_max, ind] = max(abs(c_ij));
        
        refj_rev = reverse(Y(:,j));
        c_ij2    = (circcorr(refi, refj_rev));
        [c_rev_mac, ind_rev] = max(abs(c_ij2));
        
        if c_max>c_rev_mac
            ang =(ind-1)*base_ang;
            relative_actions(indi,indj) = mat(ang);
        else
            ang = (ind_rev-1)*base_ang;
            relative_actions(indi,indj) = mat(ang)*[1 0;0 -1];
        end
    end
end

% apply synchronization
g_array = sync_D2n(relative_actions, L);

% aligning and averaging
est = zeros(L,1);
for j=1:N
    aligned_y = apply_group_action_D2n(g_array(:,:,j)', Y(:,j));
    est = est + aligned_y;
end
est = est/N;

end
