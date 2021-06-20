% run_sync_comparison_server
%
% This script compares methods for MRA over dihedral group D_2L
% The methods are:
% 1 - cc+sync
% 2 - MoM based on LS 
% 3 - EM
% We run two cases and check error rates and running time
%
% -----------------

clear
close all

% test parameters  % => server settings
to_save = 0;         % =>1
parallel_pool = 0;   % =>1

max_repeated_trials = 1; % =>50;
num_of_noise_vals   = 5; % =>10;

basic_sigma_vals = logspace(-1.5, 0, num_of_noise_vals); 
N_values = [100, 1000]; 

% signal length, also the group parameter
L = 15;

% converting the group elements to mat form
base_ang = 2*pi/L;
mat = @(x) [cos(x),-sin(x); sin(x), cos(x)];

% the least squares is based on multiply initial guesses
num_initial_guess_LS = 10;

% EM parameters
tol      = 1e-4;
max_iter = 400;

% parallel run
if parallel_pool
    parallel_nodes = 10;    % <=== note this one and your machine's capabilities
    if isempty(gcp('nocreate'))
        parpool(parallel_nodes, 'IdleTimeout', 240);
    end
end

%%
if to_save
    folder_name = ['MRA_D2n_comparison_',datestr(now,'mmmm_dd_yy_HHMM')];
    mkdir(folder_name)
    cd(folder_name)
end

% initialization
total_snr_level = zeros(length(basic_sigma_vals),length(N_values));

total_err_Sync = zeros(length(basic_sigma_vals),length(N_values));
total_err_LS     = zeros(length(basic_sigma_vals),length(N_values));
total_err_EM    = zeros(length(basic_sigma_vals),length(N_values));
total_time_sync = zeros(length(basic_sigma_vals),length(N_values));
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
    
    err_Sync = zeros(length(basic_sigma_vals),max_repeated_trials);
    err_LS     = zeros(length(basic_sigma_vals),max_repeated_trials);
    err_EM    = zeros(length(basic_sigma_vals),max_repeated_trials);
    
    time_Sync = zeros(length(basic_sigma_vals),max_repeated_trials);
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
        
        % ground truth group elements in matrix form
        groundtruth = zeros(2,2,N);
        for i=1:N
            current_shitf = shifts(i)+1;  % normalize the shifts to group index
            if current_shitf<L
                ang =(current_shitf-1)*base_ang;
                groundtruth(:,:,i) = mat(ang);
            else
                ang = (current_shitf-L-1)*base_ang;
                groundtruth(:,:,i) = mat(ang)*[1 0;0 -1];
            end
        end
        % ---- debuging ----
        % the groundtruth ratios
        %     GT = zeros(2*N);
        %     for i=1:N
        %         indi = (2*(i-1)+1):(2*i);
        %         for j=(i+1):N
        %             indj = (2*(j-1)+1):(2*j);
        %             GT(indi,indj) = groundtruth(:,:,i)*(groundtruth(:,:,j)');
        %         end
        %     end
        % ---- debuging ----
        
        % Second LOOP: changing the noise level
        for noise_level = 1:length(sigma_vals)
            
            % the current data
            added_noise = sigma_vals(noise_level)*noise_term;
            yj          = X + added_noise;
            snr_level(noise_level, trialNum) = x_square/(norm(added_noise(:,1))^2);
            
            % Applying sync
            tic
            est = MRA_D2n_sync(yj);
            time_Sync(noise_level, trialNum) = toc;
            
            % the relative error: sync
            err_Sync(noise_level, trialNum) = relative_error_D2n(est, x);
            
            % Applying LS
%             if snr_level(noise_level, trialNum)<1.5
%                 nig = num_initial_guess_LS;
%             else
%                 nig = 4;
%             end
            
            warning('off','all')
            nig = num_initial_guess_LS;
            tic
            [x_est_ls, est_dist_ls] = MRA_dihedral_LS_multi_guess(yj, sigma_vals(noise_level), nig);
            time_LS(noise_level, trialNum) = toc;
            warning('on','all')
            
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
    
    total_err_Sync(:,n_val)= mean(err_Sync,2);
    total_err_LS(:,n_val)  = mean(err_LS,2);
    total_err_EM(:,n_val)  = mean(err_EM,2);
    
    total_time_sync(:,n_val)  = mean(time_Sync,2);
    total_time_LS(:,n_val)  = mean(time_LS,2);
    total_time_EM(:,n_val)  = mean(time_EM,2);
    %% infer the result -- a comparison
    
%     figure;   % relative error - loglog MEDIAN
%     set(0,'defaulttextInterpreter','latex')
%     loglog(total_snr_level(:,n_val), median(err_Sync,2),'-.b','LineWidth',3.5);
%     hold on;
%     loglog(total_snr_level(:,n_val), median(err_LS,2),'r','LineWidth',3.3);
%     loglog(total_snr_level(:,n_val), median(err_EM,2),':k','LineWidth',3.3);
%     loglog(total_snr_level(:,n_val), ones(size(total_snr_level(:,n_val))),'--r','LineWidth',1.5);
%     xlabel('SNR')
%     ep = 0.01;
%     xlim([min(total_snr_level(:,n_val))-ep,max(total_snr_level(:,n_val))+ep])
%     ylabel('Relative error')
%     leg1 = legend('Sync','LS MoM','EM','location','best'); %best');
%     set(leg1,'Interpreter','latex');
%     set(leg1,'FontSize',20);
%     set(gca,'FontSize',18)
%    % title(['Number of elements: ', num2str(N)]);
%     if to_save
%         name_it = ['median_error_N_',num2str(N_values(n_val))];
%         saveas(gcf, name_it ,'fig');
%         saveas(gcf, name_it ,'jpg');
%          print(gcf, name_it ,'-dpng','-r300');
%         print('-depsc2',name_it);
%         print('-depsc2',name_it);
%     end
    
    figure;   % relative error - loglog
    set(0,'defaulttextInterpreter','latex')
    loglog(total_snr_level(:,n_val), total_err_Sync(:,n_val),'-.b','LineWidth',3.5);
    hold on;
    loglog(total_snr_level(:,n_val), total_err_LS(:,n_val),'r','LineWidth',3.3);
    loglog(total_snr_level(:,n_val), total_err_EM(:,n_val),':k','LineWidth',3.3);
   % loglog(total_snr_level(:,n_val), ones(size(total_snr_level(:,n_val))),'--r','LineWidth',1.5);
    xlabel('SNR')
    ep = 0.01;
    xlim([min(total_snr_level(:,n_val))-ep,max(total_snr_level(:,n_val))+ep])
    ylabel('Relative error')
    leg1 = legend('Synchronization','Method of moments','EM','location','NorthEast'); %best');
    set(leg1,'Interpreter','latex');
    set(leg1,'FontSize',20);
    set(gca,'FontSize',18)
   % title(['Number of elements: ', num2str(N)]);
    if to_save
        name_it = ['error_N_',num2str(N_values(n_val))];
        saveas(gcf, name_it ,'fig');
        saveas(gcf, name_it ,'jpg');
         print(gcf, name_it ,'-dpng','-r300');
        print('-depsc2',name_it);
        print('-depsc2',name_it);
    end
   
    figure;   % runtime
    set(0,'defaulttextInterpreter','latex')
    semilogx(total_snr_level(:,n_val), total_time_sync(:,n_val),'-.b','LineWidth',3.5);
    hold on;
    semilogx(total_snr_level(:,n_val), total_time_LS(:,n_val),'r','LineWidth',3.3);
    semilogx(total_snr_level(:,n_val), total_time_EM(:,n_val),':k','LineWidth',3.3);
    xlabel('SNR')
    ep = 0.01;
    xlim([min(total_snr_level(:,n_val))-ep,max(total_snr_level(:,n_val))+ep])
    ylabel('Time (sec.)')
    leg1 = legend('Synchronization','Method of moments','EM','location','NorthEast'); %best');
    set(leg1,'Interpreter','latex');
    set(leg1,'FontSize',20);
    set(gca,'FontSize',18)
   % title(['Number of elements: ', num2str(N)]);
    if to_save
        name_it = ['time_N_',num2str(N_values(n_val))];
        saveas(gcf, name_it ,'fig');
        saveas(gcf, name_it ,'jpg');
        print(gcf, name_it ,'-dpng','-r300');
        print('-depsc2',name_it);
        print('-depsc2',name_it);
    end
end

if to_save
    save('run_sync_D2n_data');
    cd '../'
end