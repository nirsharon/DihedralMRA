% script name: "gradient_check" using "ManOpt"

% ---------------------------------------------
% NOTE: one needs MANOPT (https://www.manopt.org/) for this test
% ---------------------------------------------

clear
clc; 
close all; 
tic

%% parameters
L = 10;   % signal length
N = 10^5;
sigma = 0; %.01;

%% data
x  = randn(L, 1);   % underlying signal

rho = rand(2*L, 1); % rand distribution
rho = rho/sum(rho); % raw stochastic

[X, shifts] = generate_observations(x, N, sigma, rho); 
mu = mean(X,2);
M  = (1/N)*X*(X')-sigma^2*eye(L);
lambda = 1; 

%% Create the problem structure.
%manifold = euclideancomplexfactory(L,1);
manifold = euclideanfactory(3*L,1);

problem.M = manifold;

problem.costgrad  = @(x) LS_cost_grad(x(1:L), x(1+L:end), mu, M, lambda, sigma);  

%problem.hess = TBA
%options.HessianFcn = [];

%% Numerically check gradient consistency
checkgradient(problem);
%checkhessian(problem);