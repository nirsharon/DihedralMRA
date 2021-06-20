function [ind] = mat2ind(n, mat)
% convert representation in D_2n

detflag = 0;
if det(mat)<0
    mat= mat*[1 0;0 -1];
    detflag = 1;
end
ang = mod(atan2(mat(2,1),mat(1,1))+2*pi,2*pi);
ind  = round(ang*n/(2*pi))+1;
%ang_num  = floor(unknown_ang*n/(2*pi));

if detflag == 1
    ind = ind + n;
end   

end

