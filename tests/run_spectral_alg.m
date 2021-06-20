% script name: "run_spectral_alg"
% Running the spectral algorithm
%%%%%%%%%%% EVEN CASE %%%%%%%%

clear
clc; 
close all; 

%% parameters
L = 10;   % signal length
N = 10^6; % number of observations
sigma = 0.01; %.01;

%% data
true_x   = randn(L, 1);    % underlying signal

true_rho = rand(2*L, 1);   % rand distribution
true_rho = true_rho/sum(true_rho); % raw stochastic
p = true_rho(1:L); q = true_rho(L+1:end);  % notational 

% imposing the assumption on the distribution
if mod(L,2)==0    % even length
    split = rand;
    q(1:2:end) = rand;
    q(2:2:end) = rand;
    q = q/(sum(q)/(1-sum(p)));
    q0 = q(1); 
    q1 = q(2); 
    
    p = (p + reverse(p))/2;
end
true_rho(L+1:end) = q;
true_rho(1:L) = p;

[X, shifts] = generate_observations(true_x, N, sigma, true_rho); 

%% fundamentals

circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';
Cx = circ(true_x);

F  = 1/(sqrt(L))*fft(eye(L));
D_Fx = diag(F*true_x);

M2 = Cx*diag(p)*Cx.' + Cx.'*diag(q)*Cx;
M2_fourier = F*M2*F';


%% run spectral algorithm
[est_x, ~] = spectral_alg(X, sigma, true_x, true_rho);

%% plotting 
x_est_al = align_to_reference(est_x, true_x);
err = norm(x_est_al - true_x)/norm(true_x);
fprintf('Relative error = %.4g\n', err);

figure; 
%subplot(121);
hold on; 
stem(true_x); 
stem(x_est_al);
legend({'$x$', '$\hat{x}$ (est)'},'Interpreter','latex');

% subplot(122);
% stem(true_rho); 
% stem(est_dist);
% legend('\rho', '\hat{\rho}'); 
