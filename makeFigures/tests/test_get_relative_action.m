% script name: "test_get_relative_action"

close all;
clear

L = 15;
N = 5;

% the signal 
x = randn(L,1);

noise_lev = [0, .5];
n = length(noise_lev);

shift = mod(randi(L),L)+1;
y     = circshift(x, shift);

noise_term = randn(L,1);

%%
figure;
plot(x,'LineWidth',3.2);
hold on;
plot(y,'-.','LineWidth',3.5);
leg1 = legend('Signal','Shifted copy'); %,'Location','Best');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);
%%
figure;
for j=1:n
    base_ind = (j-1)*3;
    subplot(n,3,base_ind+1);
    xj = x + noise_term*noise_lev(j);
    plot(xj,'r','LineWidth',3); 
    ylim([-max(abs(xj)),max(abs(xj))])
    
    subplot(n,3,base_ind+2); 
    yj = y + noise_term*noise_lev(j);
	plot(yj,'r','LineWidth',3);   
    ylim([-max(abs(yj)),max(abs(yj))])

    subplot(n,3,base_ind+3);
    c = Reflection(circcorr(xj, yj));
    cc = abs(c)/max(abs(c));
    stem(cc ,'filled');
    hold on;
    stem(shift+1,1,'r',':diamondr')
    [~, indm] = max(abs(c));
    if indm~=(shift+1)
        ths = ones(L,1)*c(shift+1)/max(abs(c));
        plot(ths,':r');
     %   legend('C-C','The shift','Threshold','Location','northeastoutside') 
    else
      %  legend('C-C','The shift','Location','best') 
    end
    if j==1
          legend('C-C','The shift','Location','best') 
    end
end

%%
figure;
plot(xj,'LineWidth',3.2);
hold on;
plot(yj,'-.','LineWidth',3.5);
leg1 = legend('Noisy signal','Noisy shifted copy'); %,'Location','Best');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',18);

%%
snr_level = zeros(n,1);
for j=2:n
    %snr(x,noise_term*noise_lev(j))
    snr_level(j) = norm(x)^2/norm(noise_term*noise_lev(j))^2;
end
%snr_level

%%
figure;
for j=2:n
    subplot(1,n-1,j-1);
   % subplot(2,2,j-1);
    xj = x + noise_term*noise_lev(j);
    yj = y + noise_term*noise_lev(j);
    c = Reflection(circcorr(xj, yj));
    cc = abs(c)/max(abs(c));
    stem(cc ,'filled');
    hold on;
    stem(shift+1,1,'r',':diamondr')
    title(['SNR = ', num2str(snr_level(j),'%1.1f')],'interpreter','latex')
    [~, indm] = max(abs(c));
    if indm~=(shift+1)
        ths = abs(ones(L,1)*c(shift+1)/max(abs(c)));
        plot(ths,':r','LineWidth',2.5);
    end
    if j==2
          legend('C-C','The shift','Location','best') 
    end 
end
