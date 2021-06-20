% test_sync_D2n
%
% A simple test for the method based on synchronization

clear;

n = 8;
N = 10;

[GT, H] = sync_data_for_testing(n, N);

% clean comparison
[est] = sync_D2n(H, n);
[err, global_g] = sync_error_D2n(est, GT, n);
disp(['clean measurements, absolute error is: ',num2str(err)])

% dummy noisy comparison
noise_mat = rand(size(H));
noise_mat = (noise_mat + noise_mat')/(5*norm(noise_mat)) + H;

[est_noise] = sync_D2n(noise_mat, n);
[err2, global_g] = sync_error_D2n(est_noise, GT, n);
disp(['noisy measurements, absolute error is: ',num2str(err2)])
