function [m] = ind2mat(n, ind)
%
base_ang = 2*pi/n;
mat = @(x) [cos(x),-sin(x); sin(x), cos(x)];

if ind<=n
    ang =(ind-1)*base_ang;
    m = mat(ang);
else
    j = ind-n;
    ang = (j-1)*base_ang;
    m = mat(ang)*[1 0;0 -1];
end

end