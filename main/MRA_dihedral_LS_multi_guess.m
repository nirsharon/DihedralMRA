function [x_est, est_dist] = MRA_dihedral_LS_multi_guess(X, sigma, num_of_guesses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Dihedral MRA: Least Squares solver 
%
% Solving MRA with LS fitting over the dihedral group
% We are using both first and second moments with a weight on the latter 
% which depends on the noise level.
% The "initial_guess" variable is optional. The defualt is a random guess.
%
% INPUT
% num_of_guesses
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<3
    num_of_guesses = 5;
end

% basic parameters
[L,N] = size(X);
if ~exist('sigma', 'var') || isempty(sigma)
    sigma = std(sum(X, 1))/sqrt(L);      % sample sigma
end

% the moments
mu = mean(X,2);
M = (1/N)*X*(X'); % we debias inside the cost function 

% the cost function weight 
lambda = 1/(L*(1+sigma^2));  % 1/(L*(1+3*sigma^2));

% cost function
fun = @(x) LS_cost_grad(x(1:L), x(1+L:end), mu, M, lambda, sigma);

% constraints
A = blkdiag(zeros(L), -eye(2*L));
b = zeros(3*L,1);
Aeq = [zeros(1,L),ones(1,2*L)];
beq = 1;

% solver options
options = optimoptions(@fmincon,'Algorithm','interior-point');
options.MaxIterations          = 2000;
options.MaxFunctionEvaluations = 10000;
options.FunctionTolerance = 1e-25;
options.ConstraintTolerance =  1e-8;
options.StepTolerance =  1e-18; 
options.Display = 'off';
options.SpecifyObjectiveGradient = true(1);
options.OptimalityTolerance = 1e-20; 
options.FiniteDifferenceStepSize = 1e-25;

% main loop
est = zeros(3*L, num_of_guesses);
fval = zeros(num_of_guesses, 1);

for trial_num = 1:num_of_guesses
    
    % current initial guess
    init_dist = rand(2*L,1); 
    init_dist = init_dist/sum(init_dist);
    init_x = rand(L,1); 
    init_x = init_x*(mean(mu)/mean(init_x)); %set the scale
    x0 = [init_x ;init_dist];
    
    % solving
    [est(:,trial_num),  fval(trial_num), exitflag2, output2]= fmincon(fun,x0,A,b,Aeq,beq,[],[],[],options);
end

% concluding
[~, ind] = min(fval);
x_est      = est(1:L, ind);
est_dist   = est((L+1):3*L, ind); 
x_est      = x_est(:);
est_dist   = est_dist(:);

end