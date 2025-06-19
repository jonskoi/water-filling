clear all;
clc;

%% parameters
D = 1000;    % number of experiments
N = 4;    % subcarriers
E = 100;   % total energy => Variable

%% water-filling algorithm loop
for i = 1:D
    g_k = abs(randn(1,N) + 1j*randn(1,N)).^2 ;
    g_k_sort = sort(g_k,'descend');
    g_k_sort_inv = 1./g_k_sort;                     % sort g_k and take inverse of it

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

    data_rate(i) = sum(data_rate_sc);               % total data rate for single experiment
    data_rate_sc = 0;

end

avg_data_rate = 1/D*sum(data_rate)                 % average of data rate = Expected data rate

%% fig 1: energy allocation example from 1000th experiment
figure;
subcarrier_index = 1:N; 
mu_values = mu * ones(1, N); 
eps_plot = eps;
g_k_sort_inv_plot = 1./g_k_sort;                  % inverse of g_k before elimination
bar_data = [g_k_sort_inv_plot; eps_plot]';          % for bar graph

hold on;
bar(subcarrier_index, bar_data, 'stacked'); 
plot(subcarrier_index, mu_values, 'LineWidth', 2, 'Color', 'black');

xlabel('Subcarrier Index');
xticks(1:N);
legend('1/g_k', '\epsilon_k', '\mu', 'Location', 'NorthEast');
title(['\epsilon = ', num2str(E)]);
grid on;
hold off;