function [err, global_element] = sync_error_D2n(A, B, n)
% 
% calculate the D2n global ambiguity to minimize the error between A and B
%
% April 21

N = size(A,3);

% if det(A(:,:,1))~=det(B(:,:,1))
%     for j=1:N
%         A(:,:,j) = A(:,:,j)*[1 0;0 -1]; 
%     end
% end

Q = A(:,:,1)'*B(:,:,1);
    
diff = zeros(size(A));
for j=1:N
    diff(:,:,j)  = A(:,:,j)*Q*B(:,:,j)';
end

P = sum(diff,3);
[ u, ~, v] = svd(P);
ortho = u*v';

% rounding -- to D_2n
[global_element] = D2n_rounding(ortho, n);
% detflag = 0;
% if det(global_element)<0
%     global_element = [1 0;0 -1]*global_element;
%     detflag = 1;
% end
% base_ang = mod(atan2(global_element(2,1),global_element(1,1))+2*pi,2*pi);
% ang = floor(base_ang*n/(2*pi));
% 
% global_element = [cos(ang),-sin(ang); sin(ang), cos(ang)];
% if detflag == 1
%     global_element = [1 0;0 -1]*global_element;      
% end   

global_element = global_element*Q;

err = 0;
for j=1:N
    err = err + norm(A(:,:,j)*global_element-B(:,:,j),'fro')^2;  % squared 
end

% averaged squared error
err = err/N;

end

