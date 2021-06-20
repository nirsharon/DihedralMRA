%script name: "test_D2n_rounding"

clear;
n = 5;

base_ang = 2*pi/n;
mat = @(x) [cos(x),-sin(x); sin(x), cos(x)];


% test 1: returning the same element
j = n-2;
ang =(j-1)*base_ang;
g = mat(ang);
norm(g - D2n_rounding(g, n))

% test 2: go in-between to see the rounding
j = n-2+.51;
ang =(j-1)*base_ang;
g = mat(ang);
norm(g - D2n_rounding(g, n))