% script name: "addPath"
% Setting the path for the DihedralMRA package
%

fprintf('Adding path: \n');

addpath(genpath('main'));
addpath(genpath('auxiliaries'));
addpath(genpath('tests'));

x = input('Do you interest in reproducing figures from the paper? \n press ''y'' for yes. Otherwise, press any key \n','s');
if strcmp(x, 'y')
    addpath(genpath('makeFigures'));
    fprintf(' "/makeFigures" is added \n');
end

fprintf('Done! \n');

