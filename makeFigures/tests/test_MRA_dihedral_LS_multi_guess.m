% script name: "test_MRA_dihedral_LS_multi_guess"
%
% testing the multi-guess version of the MoM LS

clear
clc; 
close all; 

tic

%% parameters
L = 11;   % signal length
N = 10^6; % number of observations
sigma = 0.001; %.01;

%% data
true_x   = randn(L, 1);    % underlying signal
true_rho = rand(2*L, 1);   % rand distribution
true_rho = true_rho/sum(true_rho); % raw stochastic
[X, shifts] = generate_observations(true_x, N, sigma, true_rho); 

%% comparison run
[x_est, est_dist] = MRA_dihedral_LS_multi_guess(X, sigma);
[est_x_oneGuess, est_dist_oneGuess] = MRA_dihedral_LS(X, sigma);


%% plotting 
x_est_al  = align_to_reference(x_est, true_x);
x_est_al1 = align_to_reference(est_x_oneGuess, true_x);

err = norm(x_est_al - true_x)/norm(true_x);
err1 = norm(x_est_al1 - true_x)/norm(true_x);

fprintf('Relative error multi guesses = %.4g\n', err);
fprintf('Relative error one guess = %.4g\n', err1);


figure; 
hold on; 
stem(true_x); 
stem(x_est_al);
stem(x_est_al1);
legend({'Original', 'Multi guesses', 'One guess'},'Interpreter','latex');
set(gca,'FontSize',18)

figure; 
hold on; 
stem(true_rho); 
stem(est_dist);
legend({'True dist.', 'Estimation'},'Interpreter','latex');
set(gca,'FontSize',18)

toc()

