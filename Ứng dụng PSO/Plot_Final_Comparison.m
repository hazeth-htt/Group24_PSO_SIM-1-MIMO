% Plot_Ultimate_Comparison.m
clc; clearvars; close all;

K_values = [1, 2, 5, 10];
Colors = {'r', 'g', 'b', 'm'};

% HÌNH 1: CAPACITY
figure(1); hold on; box on; grid on; set(gcf, 'Position', [100, 100, 900, 600]);
for i = 1:length(K_values)
    K = K_values(i);
    
    % 1. Benchmark (Nét đứt mỏng)
    if exist(sprintf('Results/Benchmark_Capacity_K%d.mat', K), 'file')
        load(sprintf('Results/Benchmark_Capacity_K%d.mat', K));
        plot(1:10, Capacity_average, [Colors{i} '--'], 'LineWidth', 1.5, 'DisplayName', sprintf('Benchmark K=%d', K));
    end
    
    % 2. PSO Thường (Nét liền + Tròn)
    if exist(sprintf('Results/PSO_Capacity_K%d.mat', K), 'file')
        load(sprintf('Results/PSO_Capacity_K%d.mat', K));
        plot(1:10, Capacity_average, [Colors{i} '-o'], 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', sprintf('PSO K=%d', K));
    end
    
    % 3. DEPSO Biến thể (Nét liền đậm + Sao) -> "Ngôi sao sáng"
    if exist(sprintf('Results/DEPSO_Capacity_K%d.mat', K), 'file')
        load(sprintf('Results/DEPSO_Capacity_K%d.mat', K));
        plot(1:10, Capacity_average, [Colors{i} '-*'], 'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', sprintf('DEPSO (Variant) K=%d', K));
    end
end
xlabel('Layers (L)'); ylabel('Sum Rate (bits/s/Hz)'); title('Comparison: Benchmark vs PSO vs DEPSO'); legend('Location', 'best');

% HÌNH 2: NMSE
figure(2); hold on; box on; grid on; set(gcf, 'Position', [150, 150, 900, 600]);
for i = 1:length(K_values)
    K = K_values(i);
    
    % Vẽ tương tự cho NMSE
    if exist(sprintf('Results/Benchmark_NMSE_K%d.mat', K), 'file')
        load(sprintf('Results/Benchmark_NMSE_K%d.mat', K));
        plot(1:10, NMSE_average, [Colors{i} '--'], 'LineWidth', 1.5, 'DisplayName', sprintf('Benchmark K=%d', K));
    end
    if exist(sprintf('Results/PSO_NMSE_K%d.mat', K), 'file')
        load(sprintf('Results/PSO_NMSE_K%d.mat', K));
        plot(1:10, NMSE_average, [Colors{i} '-o'], 'LineWidth', 1.5, 'DisplayName', sprintf('PSO K=%d', K));
    end
    if exist(sprintf('Results/DEPSO_NMSE_K%d.mat', K), 'file')
        load(sprintf('Results/DEPSO_NMSE_K%d.mat', K));
        plot(1:10, NMSE_average, [Colors{i} '-*'], 'LineWidth', 2.5, 'DisplayName', sprintf('DEPSO K=%d', K));
    end
end
xlabel('Layers (L)'); ylabel('NMSE'); title('NMSE Comparison'); legend('Location', 'best');