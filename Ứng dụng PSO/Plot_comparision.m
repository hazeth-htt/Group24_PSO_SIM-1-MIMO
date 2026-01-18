% Plot_Comparison_Full.m
% Code vẽ so sánh đầy đủ giữa Benchmark (Gradient) và PSO (Proposed)
% Vẽ cả Capacity và NMSE, sau đó lưu ảnh tự động.

clc; clearvars; close all;

%% CẤU HÌNH VẼ
K_values = [1, 2, 5, 10];
Colors = {'r', 'g', 'b', 'm'}; % Đỏ, Xanh lá, Xanh dương, Tím
LineStyles_Bench = '--';       % Benchmark: Nét đứt
LineStyles_PSO = '-';          % PSO: Nét liền
Markers_PSO = {'o', 's', '^', 'd'}; % Các loại dấu cho PSO dễ nhìn

% Kiểm tra thư mục kết quả
if ~exist('Results', 'dir')
    error('Khong tim thay thu muc Results! Hay chay 2 file test_SIM truoc.');
end

%% --- HÌNH 1: SO SÁNH DUNG LƯỢNG (CAPACITY) ---
figure(1); hold on; box on; grid on;
set(gcf, 'Position', [100, 100, 800, 600]); % Chỉnh kích thước cửa sổ to đẹp

% 1. Vẽ Benchmark
for i = 1:length(K_values)
    K = K_values(i);
    FileName = sprintf('Results/Benchmark_Capacity_K%d.mat', K);
    if exist(FileName, 'file')
        load(FileName);
        plot(1:10, Capacity_average, [Colors{i} LineStyles_Bench], ...
            'LineWidth', 2, 'DisplayName', sprintf('Benchmark (K=%d)', K));
    end
end

% 2. Vẽ PSO
for i = 1:length(K_values)
    K = K_values(i);
    FileName = sprintf('Results/PSO_Capacity_K%d.mat', K);
    if exist(FileName, 'file')
        load(FileName);
        plot(1:10, Capacity_average, [Colors{i} LineStyles_PSO Markers_PSO{i}], ...
            'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', Colors{i}, ...
            'DisplayName', sprintf('Proposed PSO (K=%d)', K));
    end
end

% Trang trí Hình 1
xlabel('Number of Transmit Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Sum Rate (bits/s/Hz)', 'FontSize', 12, 'FontWeight', 'bold');
title('Comparison of Sum-Rate Capacity', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
saveas(gcf, 'Results/Comparison_Capacity.png'); % Lưu ảnh

%% --- HÌNH 2: SO SÁNH SAI SỐ (NMSE) ---
figure(2); hold on; box on; grid on;
set(gcf, 'Position', [150, 150, 800, 600]);

% 1. Vẽ Benchmark
for i = 1:length(K_values)
    K = K_values(i);
    FileName = sprintf('Results/Benchmark_NMSE_K%d.mat', K);
    if exist(FileName, 'file')
        load(FileName);
        plot(1:10, NMSE_average, [Colors{i} LineStyles_Bench], ...
            'LineWidth', 2, 'DisplayName', sprintf('Benchmark (K=%d)', K));
    end
end

% 2. Vẽ PSO
for i = 1:length(K_values)
    K = K_values(i);
    FileName = sprintf('Results/PSO_NMSE_K%d.mat', K);
    if exist(FileName, 'file')
        load(FileName);
        plot(1:10, NMSE_average, [Colors{i} LineStyles_PSO Markers_PSO{i}], ...
            'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', Colors{i}, ...
            'DisplayName', sprintf('Proposed PSO (K=%d)', K));
    end
end

% Trang trí Hình 2
xlabel('Number of Transmit Layers (L)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('NMSE (Normalized Mean Square Error)', 'FontSize', 12, 'FontWeight', 'bold');
title('Comparison of Channel Fitting NMSE', 'FontSize', 14);
legend('Location', 'best', 'FontSize', 10);
saveas(gcf, 'Results/Comparison_NMSE.png'); % Lưu ảnh

fprintf('Da ve va luu xong 2 hinh vao thu muc Results!\n');