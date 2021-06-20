% check LS

clear
close all
to_save = 0;

% signal length, also the group parameter
L = 15;

% converting the group elements to mat form
base_ang = 2*pi/L;
mat = @(x) [cos(x),-sin(x); sin(x), cos(x)];

% test parameters
max_repeated_trials = 1;
basic_sigma_vals = logspace(-3.5, .7 ,10); 
N_values = 10^4;

% the least squares is based on multiply initial guesses
num_initial_guess_LS = 10;

% EM parameters
tol      = 1e-3;
max_iter = 200;

circ = @(v) toeplitz([v(1); v(end:-1:2)], v)';

% parallel run
parallel_pool = 1;
if parallel_pool
    parallel_nodes = 10;
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

%%
if to_save
    folder_name = ['MRA_D2n_LargeN_',datestr(now,'mmmm_dd_yy_HHMM')];
    mkdir(folder_name)
    cd(folder_name)
end

% initialization
total_snr_level = zeros(length(basic_sigma_vals),length(N_values));
total_err_LS     = zeros(length(basic_sigma_vals),length(N_values));
total_err_EM    = zeros(length(basic_sigma_vals),length(N_values));
total_time_LS     = zeros(length(basic_sigma_vals),length(N_values));
total_time_EM    = zeros(length(basic_sigma_vals),length(N_values));

for n_val = 1:length(N_values)
    N = N_values(n_val);
    
    % defining the signal
    x = randn(L,1);
    x_square = norm(x)^2;
    
    % defining the distribution (for MoM LS)
    p = rand(2*L, 1);
    p = p/sum(p);
    
    % initialization
    snr_level = zeros(length(basic_sigma_vals),max_repeated_trials);
    
    err_LS     = zeros(length(basic_sigma_vals),max_repeated_trials);
    err_EM    = zeros(length(basic_sigma_vals),max_repeated_trials);
    
    time_LS     = zeros(length(basic_sigma_vals),max_repeated_trials);
    time_EM    = zeros(length(basic_sigma_vals),max_repeated_trials);
    
    % Repeated trial LOOP
    for trialNum = 1:max_repeated_trials      
        
        % defining a unified noise term for eac trial
        noise_term = randn(L,N);  %rand(L,N);
        noise_fact   = 1.2*sqrt(x_square/norm(noise_term(:,1))^2);
        sigma_vals = basic_sigma_vals*noise_fact;
        
        % make MRA data: fixed for each trial (no noise yet)
        [X, shifts] = generate_observations(x, N, 0, p);
        
        % Second LOOP: changing the noise level
        for noise_level = 1:length(sigma_vals)
            
            % the current data
            added_noise = sigma_vals(noise_level)*noise_term;
            yj          = X + added_noise;
            snr_level(noise_level, trialNum) = x_square/(norm(added_noise(:,1))^2);
            
          
            % Applying LS    
            warning('off','all')
            nig = num_initial_guess_LS;
            tic
            [x_est_ls, est_dist_ls] = MRA_dihedral_LS_multi_guess(yj, sigma_vals(noise_level), nig);
            time_LS(noise_level, trialNum) = toc;
            warning('on','all')
            
            Cx = circ(x);
            pp = p(1:L);
            qq = p(L+1:end);
            M1 = Cx*pp  + Cx'*qq; 
            M2 = Cx*diag(pp)*Cx.' + Cx.'*diag(qq)*Cx;
            mu1 = mean(yj,2);
            mu2 = (1/N)*yj*(yj')-sigma_vals(noise_level)^2*eye(L);
            norm(mu1(:)-M1(:))
            norm(mu2(:)-M2(:))
            
            % the relative error: LS
            err_LS(noise_level, trialNum) = relative_error_D2n(x_est_ls, x);
            
            % initial guess for the EM
            x_init = randn(L, 1);
            p_init = rand(2*L, 1);
            p_init = p_init/sum(p_init);
            
            % Applying EM
            tic
            [x_est_em, p_est_em, likelihood] = EM(yj, sigma_vals(noise_level), x_init, p_init, max_iter, tol);
            time_EM(noise_level, trialNum) = toc;
            
            % the relative error: EM
            err_EM(noise_level, trialNum) = relative_error_D2n(x_est_em, x);
        end
        % trialNum   % print the trial number for control
    end
    
    %% summary
    total_snr_level(:,n_val) = mean(snr_level,2);
    total_err_LS(:,n_val)  = mean(err_LS,2);
    total_err_EM(:,n_val)  = mean(err_EM,2);
    
    total_time_LS(:,n_val)  = mean(time_LS,2);
    total_time_EM(:,n_val)  = mean(time_EM,2);
    %% infer the result -- a comparison     
    
    figure;   % runtime
    set(0,'defaulttextInterpreter','latex')
    semilogx(total_snr_level(:,n_val), total_time_LS(:,n_val),'r','LineWidth',3.3);
    hold on;
    semilogx(total_snr_level(:,n_val), total_time_EM(:,n_val),':k','LineWidth',3.3);
    xlabel('SNR')
    ep = 0.01;
    xlim([min(total_snr_level(:,n_val))-ep,max(total_snr_level(:,n_val))+ep])
    ylabel('Time (sec.)')
    leg1 = legend('LS MoM','EM','location','NorthEast'); 
    set(leg1,'Interpreter','latex');
    set(leg1,'FontSize',20);
    set(gca,'FontSize',18)
    %title(['Number of elements: ', num2str(N)]);
    if to_save
        name_it = ['comp_time_D2n_N_',num2str(n_val)];
        saveas(gcf, name_it ,'fig');
        saveas(gcf, name_it,'jpg');
        print('-depsc2',name_it);
        print('-depsc2',name_it);
    end
    
    figure;   % relative error - loglog
    set(0,'defaulttextInterpreter','latex')
    loglog(total_snr_level(:,n_val), total_err_LS(:,n_val),'r','LineWidth',3.3);
    hold on;
    loglog(total_snr_level(:,n_val), total_err_EM(:,n_val),':k','LineWidth',3.3);
    loglog(total_snr_level(:,n_val), ones(size(total_snr_level(:,n_val))),'--r','LineWidth',1.5);
    xlabel('SNR')
    ep = 0.01;
    xlim([min(total_snr_level(:,n_val))-ep,max(total_snr_level(:,n_val))+ep])
    ylabel('Relative error')
    leg1 = legend('LS MoM','EM','location','NorthEast'); %best');
    set(leg1,'Interpreter','latex');
    set(leg1,'FontSize',20);
    set(gca,'FontSize',18)
   % title(['Number of elements: ', num2str(N)]);
    if to_save
        name_it = ['compare_D2n_N_',num2str(n_val)];
        saveas(gcf, name_it ,'fig');
        saveas(gcf, name_it,'jpg');
        print('-depsc2',name_it);
        print('-depsc2',name_it);
    end
end

if to_save
    save('run_compare_D2n_data');
    cd '../'
end
