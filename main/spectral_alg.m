function [est_x, est_dist] = spectral_alg(X, sigma, true_x, true_rho)
%============== Dihedral Version ============== 
% Solving MRA by a direct spectral solver. 
% ===> No distribution part yet <======
%
% Last two inputs are for debugging
%
% NS, August 20.

[L,N] = size(X);

if ~exist('sigma', 'var') || isempty(sigma)
    sigma = std(sum(X, 1))/sqrt(L);      % sample sigma
end

circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';

% estimating the two (sample) moments
mu = mean(X,2);
M = (1/N)*X*(X')-sigma^2*eye(L);

% the power spectrum (modulus of DFT of x)
P_est = (1/N)*sum(abs(fft(X)).^2,2)-L*sigma^2;
if ~isequal(1*(P_est<1e-10),zeros(L,1))
    warning('Not enough samples for power spectrum positive estimation. Projection applied.');
    P_est(P_est<1e-10) = 1e-10;
end
P_est = P_est.^(.5);

% tranform the second moment matrix M and extract eigenvectors
F             = 1/(sqrt(L))*fft(eye(L));
D_normalize_x = diag(P_est.^(-1));
conj_M        = 1/sqrt(L)*D_normalize_x*F';

% tranform the second moment matrix M for extracting eigenvectors
M_fourier = F*M*F';
vv        = zeros(L,1); 
vv(1+L/2) = 1;
index_mat = circ(vv);
ind_mat2  = ones(L)-index_mat;
M_tilde              = zeros(L);
M_tilde(ind_mat2==1) = M_fourier(ind_mat2==1);
M_ext     = conj_M'*M_tilde*conj_M;

% % ----- debuging
% Cx = circ(true_x);
% norm_mat = diag(abs(F*true_x).^(-1));
% conj_M2 = 1/sqrt(L)*norm_mat*F';
% Cx_ortho = (1/sqrt(L)*F*norm_mat*F')*Cx;
% M_ext = conj_M2'*M2_tilde*conj_M2;
% %V     = Cx_ortho'; % columns are eigenvectors

try
 [V, D] = eig(M_ext);
catch
  M_ext = (M_ext+M_ext')/2;  % impose symmetry if needed
  [V, D] = eig(M_ext);
end

% choosing distinct entry if exists 
raw_rho = abs(diag(D));
tol = 1e-4;
[ii,jj,kk] = uniquetol(raw_rho, tol);
f = histc(kk,1:numel(jj));         % frequencies corresponding to ii
idx = f<2;                         % "single" frequencies
idx1 = ismember(raw_rho,ii(idx));  % ind of the "singles"
if (1*idx1)==zeros(size(raw_rho))
    warning('spectral_method: No distinct distribution entry. Result may be inaccurate');
    % choosing eigenvector that leads to least imaginary solution
    Z = fft(V);
    [~, ind] = min(sum(imag(Z).^2));
else 
    % choosing eigenvector that leads to least imaginary solution  
    % ------- other alternative: a cost function like ----------
    % \[ norm(M_hat-circ(xhat)*diag(est_dist)*circ(xhat)') ,\]
    % where xhat = (sign(V(1,j))/sign(mean(X(:))))*V(:,j)
    % ----------------------------------------------------------
    %[~, ind_in_idx1] = min(sum(imag(fft(V(:,idx1))).^2));
   [~, ind_in_idx1] = min(sum(imag((V(:,idx1))).^2));
    val = raw_rho(idx1);    
    ind = find(~(raw_rho-val(ind_in_idx1)));
end

% estimating the signal
est_x = real(F'*(P_est.*F*V(:,ind)));

% need to restore scale;
est_x  = est_x*(mean(mu)/mean(est_x)); % scale correction

% estimating the distribution
p_and_q = (circ(est_x))\mu;
est_dist = [];

if nargin==4 % debugging mode
%     % if true values exist, output the aligned signal and rho
%     est_x = align_to_reference(est_x, true_x);
%     est_dist = align_to_reference(est_dist, true_rho);
%     figure; plot(est_x); hold on; plot(true_x);
%     norm(est_x-true_x)
%     figure; plot(est_rho); hold on; plot(true_rho);
end

end