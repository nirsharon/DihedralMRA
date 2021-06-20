function [y] = apply_group_action_D2n(g, x )
% y = g \circ x, where g \in D_2n
% g is given in matrix form
 
n = length(x);
j = mat2ind(n, g);
if j<=n
    y = circshift(x, j-1);
else
    y = circshift(reverse(x), j-1-n);
end

    

end

