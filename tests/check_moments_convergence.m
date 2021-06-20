% script name: check_moments_convergence
%
% checking the code-generated moments against the analytic formulas

clear
clc; 
close all; 
tic

%% parameters
L = 11;   % signal length
num_of_trials = 10;
N_arr = floor(logspace(2,6,num_of_trials)); % observations #, ranges from 100 to 10^5
sigma = 0.1; %.01;

%% data
x  = randn(L, 1);   % underlying signal
xr = reverse(x);    % x reflected

dist = rand(2*L, 1);   % rand distribution
dist = dist/sum(dist); % raw stochastic
p = dist(1:L); q = dist(L+1:end);  % simple notation

[X, shifts] = generate_observations(x, N_arr(end), sigma, dist); 

%% analytic moments
circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';

Cx = circ(x);
M1 = Cx*p + Cx.'*q;
M2 = Cx*diag(p)*Cx.' + Cx.'*diag(q)*Cx;


%% test convergence
err1 = zeros(num_of_trials,1);
err2 = zeros(num_of_trials,1);
for j=1:num_of_trials
    N = N_arr(j);
    currentX = X(:,1:N);
    mu1 = mean(currentX,2);
    mu2 = (1/N)*currentX*(currentX')-sigma^2*eye(L);
    err1(j) = norm(mu1(:)-M1(:));
    err2(j) = norm(mu2(:)-M2(:));
end

%% plot
figure;
loglog(N_arr, err1);
grid on;
hold on;
loglog(N_arr,exp(log(err1(1))-.5*(log(N_arr)-log(N_arr(1)))),'r' );
legend('comvergence','expected rate')
title('First moment')
set(gca,'FontSize',18)


%%
figure
loglog(N_arr, err2);
hold on;
loglog(N_arr,exp(log(err2(1))-.5*(log(N_arr)-log(N_arr(1)))),'r' );
legend('comvergence','expected rate')
grid on;
title('Second moment')
set(gca,'FontSize',18)

toc()