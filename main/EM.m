function[x, p, likelihood] = EM(X, sigma, x, p, max_iter, tol)
%
% Input:
% x -- the initial guess, if not appears then it is taken as random
% tol -- for alting criterion
%

% X contains N observations, each of length L
[L, N] = size(X);

% In practice, we iterate with the DFT of the signal x
fftx = fft(x);

% Precomputations on the observations
fftX = fft(X);
sqnormX = repmat(sum(abs(X).^2, 1), L, 1);
likelihood = zeros(max_iter, 1);

for iter = 1 : max_iter
    [fftx_new, p_new, likelihood(iter)] = EM_iteration(fftx, p, fftX, sqnormX, sigma, L, N);
    % note that abs(fft)) is invariant under the dihedral group
    if abs(abs(fftx)- abs(fftx_new)) < tol
        break;
    end
    fftx = fftx_new;
    p = p_new;
end
likelihood = likelihood(1:iter);
%fprintf('\t\tEM: %d full iterations\n', iter);
x = ifft(fftx);
end


function[fftx_new, p_new, likelihood] = EM_iteration(fftx, p, fftX, sqnormX, sigma, L, N)

C1 = ifft(bsxfun(@times, conj(fftx), fftX));
C2 = ifft(bsxfun(@times, fftx, fftX));
sqnormx = norm(fftx)^2/L;
T1 = (2*C1 - sqnormx - sqnormX)/(2*sigma^2);
T2 = (2*C2 -  sqnormx - sqnormX)/(2*sigma^2);
T = [T1; T2];
% computes the likelihood (of the previous iterate
P = repmat(p, 1, N);
likelihood = sum(log(sum(P.*exp(T),1))); 
T = bsxfun(@minus, T, max(T, [], 1));
W = exp(T);
W = bsxfun(@times, p, W);
W = bsxfun(@times, W, 1./sum(W, 1));

% the non-reflected part
fftx_new = sum(conj(fft(W(1:L,:))).*fftX, 2);
% the reflected part
fftx_new = fftx_new + conj(sum((fft(reverse(W(L+1:2*L,:)))).*(fftX), 2));
fftx_new = fftx_new/N;
p_new = mean(W,2);
end

%
% dumb way
% fftx_new2 = zeros(L,1);
% X = ifft(fftX);
% for i = 1:size(X,2)
%     vec = zeros(L, 1);
%     for l = 0:L-1
%         vec = vec + W(L+1+l,i).*reverse(circshift(X(:,i),-l),L);
%     end
%     fftx_new2 = fftx_new2 + fft(vec);
% end
%
% fftx_new = fftx_new + fftx_new2;

