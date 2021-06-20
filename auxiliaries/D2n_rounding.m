function [g_dihedral] = D2n_rounding(g, n)
%
% getting a rotation $g$, we calculate the closest D2n element
%

detflag = 0;
if det(g)<0
    g = g*[1 0;0 -1];
    detflag = 1;
end
unknown_ang = mod(atan2(g(2,1),g(1,1))+2*pi,2*pi);
ang_num  = round(unknown_ang*n/(2*pi));
%ang_num  = floor(unknown_ang*n/(2*pi));
ang = (2*pi/n)*ang_num;

g_dihedral = [cos(ang),-sin(ang); sin(ang), cos(ang)];

if detflag == 1
    g_dihedral = g_dihedral*[1 0;0 -1];      
end   

end

