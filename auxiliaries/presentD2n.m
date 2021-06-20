function [] = presentD2n(n)
%
base_ang = 2*pi/n;
mat = @(x) [cos(x),-sin(x); sin(x), cos(x)];

line1 = [];
for j=1:n
    ang =(j-1)*base_ang;
   % fprintf('%.2f %.2f\n',mat(ang))
    line1 = [line1 mat(ang)];
end
rotations = line1

line2 = [];
for j=1:n
    ang = (j-1)*base_ang;
    %fprintf('%.2f %.2f\n', [1 0;0 -1]*mat(ang))
    line2 = [line2 mat(ang)*[1 0;0 -1]];
end
reflections = line2

end

