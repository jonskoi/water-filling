clear all;
clc;

%% parameters
% D = 1000;    % number of experiments
N = 16;    % subcarriers
SNR = [-10 -5 0 5 10 15 20];  % SNR values in dB
E = 10.^(SNR / 10);           % Total energy calculated from SNR values (10^(SNR/10))

% assign one random g_k value to both algorithms
g_k = abs(randn(1,N) + 1j*randn(1,N)).^2 ;
g_k_sort = sort(g_k,'descend');
g_k_sort_inv = 1./g_k_sort; 

%% caculate data rate for both algorithms
for l = 1:length(SNR)
    [data_rate_wf(l), data_rate_eq(l),eps] = getDataRate(E(l), N, g_k, g_k_sort, g_k_sort_inv);
end

%% getDataRate function
function [data_rate_wf, data_rate_eq,eps] = getDataRate(E, N, g_k, g_k_sort, g_k_sort_inv)
    % water filling algorithm
    for n = N:-1:1
        mu = 1/n*(E + sum(g_k_sort_inv(:)));
        eps(1:n) = mu-g_k_sort_inv(1:n);            % calculate energy for each subcarriers
        
        if eps(n) < 0  
            eps(n) = 0;                             % eliminate energy allocated to the last subcarrier
            g_k_sort_inv(n) = 0;        
        else 
            break;
        end
    end

    for k = 1:n                                     % kth subcarrier    
        data_rate_sc(k) = log2(1+g_k_sort(k)*eps(k));    % calculate data rate for single subcarrier       
    end        
    data_rate_wf = sum(data_rate_sc);               % total data rate for single experiment
    data_rate_sc = 0;

    % Equal energy allocation algorithm
    eq_eps = ones(1, N) * E / N;                    % Equal energy allocation for each subcarrier
    eq_data_rate = sum(log2(1 + g_k(1:N) .* eq_eps(1:N)));  
    data_rate_eq = eq_data_rate;  
end

%% fig 2: Data rate vs SNR (dB)
figure;
hold on;
plot(SNR, data_rate_wf, '-or', 'DisplayName', 'Water-Filling');
plot(SNR, data_rate_eq, '-ob', 'DisplayName', 'Equal Energy Allocation');

% Labels and legend
xlabel('SNR (dB)');
ylabel('Data Rate (bps/Hz)');
title(['N = ', num2str(N)]);
legend('show');
grid on;
hold off;