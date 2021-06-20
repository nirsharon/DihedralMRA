function [X, shifts] = generate_observations(x, N, sigma, p)
% Given a signal x of length N, generates a matrix X of size N x M such
% that each column of X is a randomly, circularly shifted version of x with
% iid Gaussian noise of variance sigma^2 added on top.

x  = x(:);
L  = length(x);
xr = reverse(x);
X  = zeros(L, N);
rho_vec = cumsum(p);
shifts_ind = rand(N, 1);
shifts = discretize(shifts_ind, [0; rho_vec]);
shifts = shifts - 1;

for m = 1 : N
    if shifts(m)>=L %reversed signal
        X(:, m) = circshift(xr, mod(shifts(m),L));
    else
        X(:, m) = circshift(x, shifts(m));
    end
end

X = X + sigma*randn(L, N);
