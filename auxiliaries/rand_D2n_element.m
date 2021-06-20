function [g_arr] = rand_D2n_element(n, N)
% generate N uniformly random elements in D_2n

if nargin<2
    N=1;
end

g_arr    = zeros(2,2,N);
base_ang = 2*pi/n;
rand_arr = randi(n,N,1);
rand_ref = randi(2,N,1)-1;

for j=1:N
    ang = (rand_arr(j)-1)*base_ang;
    g_arr(:,:,j) = [cos(ang),-sin(ang); sin(ang), cos(ang)];
    if rand_ref(j)==0
        g_arr(:,:,j) = g_arr(:,:,j)*[1 0;0 -1];        
    end   
end


end
