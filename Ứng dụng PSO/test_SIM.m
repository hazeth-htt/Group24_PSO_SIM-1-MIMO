% test_SIM_Benchmark.m
% Code này chạy thuật toán GỐC (Gradient Descent) với K = 1, 2, 5, 10
% Đã tích hợp đầy đủ thuật toán tối ưu của tác giả.
clc; clearvars; close all;

%% CẤU HÌNH CHUNG
K_values = [1, 2, 5, 10]; % Các trường hợp cần test
MonteCarlo = 10;          % Số lần lặp (Chạy thật nên để 10-50 tùy thời gian của bạn)
Max_L = 10;               % Số lớp phát (Chạy từ 1 đến 10)

% Các thông số vật lý
Thickness = 0.05; Pt = 10^(20/10); Sigma2 = 10^(-110/10);
c = 3*10^8; f0 = 28*10^9; lambda = c/f0;
N_max = 10; PL = -20*log10(4*pi/lambda)-35*log10(250); pathloss = 10^(PL/10);
M = 100; N = 100; d_element_spacing = lambda/2; S = 4;

%% VÒNG LẶP LỚN: CHẠY TỪNG TRƯỜNG HỢP K
for k_idx = 1:length(K_values)
    K = K_values(k_idx); 
    fprintf('\n>>> DANG CHAY BENCHMARK (Gradient Descent) VOI K = %d <<<\n', K);
    
    Capacity_average = zeros(Max_L,1);
    NMSE_average = zeros(Max_L,1);
    
    % Vòng lặp số lớp phát L
    for ii = 1:Max_L
        L = ii;
        fprintf('   - Layer L=%d...', L);
        tic
        
        % --- KHỞI TẠO MA TRẬN KÊNH ---
        d_layer_spacing_transmit = Thickness/L; d_layer_spacing_receive = Thickness/K;
        W_T = zeros(M,M); Corr_T = zeros(M,M); U_R = zeros(N,N); Corr_R = zeros(N,N);
        W_T_1 = zeros(M,S); U_R_1 = zeros(S,N);
        
        % Tính toán kênh (Rút gọn)
        for mm1=1:M, for mm2=1:M, d_temp=sqrt((mod(mm1-1,10)-mod(mm2-1,10))^2+(ceil(mm1/10)-ceil(mm2/10))^2)*lambda/2; W_T(mm2,mm1)=lambda/4/pi/sqrt((Thickness/L)^2+d_temp^2)*exp(-1i*2*pi*sqrt((Thickness/L)^2+d_temp^2)/lambda); if d_temp==0, Corr_T(mm2,mm1)=1; else, Corr_T(mm2,mm1)=sin(pi*2*d_temp/lambda)/(pi*2*d_temp/lambda); end; end; end
        for nn1=1:N, for nn2=1:N, d_temp=sqrt((mod(nn1-1,10)-mod(nn2-1,10))^2+(ceil(nn1/10)-ceil(nn2/10))^2)*lambda/2; U_R(nn2,nn1)=lambda/4/pi/sqrt((Thickness/K)^2+d_temp^2)*exp(-1i*2*pi*sqrt((Thickness/K)^2+d_temp^2)/lambda); if d_temp==0, Corr_R(nn2,nn1)=1; else, Corr_R(nn2,nn1)=sin(pi*2*d_temp/lambda)/(pi*2*d_temp/lambda); end; end; end
        for mm=1:M, for nn=1:S, d_tr=sqrt((Thickness/L)^2+((mod(mm-1,10)+1-(1+10)/2)*lambda/2)^2+((ceil(mm/10)-(1+10)/2)*lambda/2-(nn-(1+S)/2)*lambda/2)^2); W_T_1(mm,nn)=lambda/4/pi/d_tr*exp(-1i*2*pi*d_tr/lambda); end; end
        for mm=1:N, for nn=1:S, d_re=sqrt((Thickness/K)^2+((mod(mm-1,10)+1-(1+10)/2)*lambda/2)^2+((ceil(mm/10)-(1+10)/2)*lambda/2-(nn-(1+S)/2)*lambda/2)^2); U_R_1(nn,mm)=lambda/4/pi/d_re*exp(-1i*2*pi*d_re/lambda); end; end

        Capacity = zeros(MonteCarlo,1); 
        NMSE = zeros(MonteCarlo,1);
        
        % Chuẩn bị biến tạm cho thuật toán Gradient
        Derivative_transmit_phase_shift = zeros(M,L);
        Derivative_receive_phase_shift = zeros(N,K);
        Temp1 = zeros(S,1); Temp2 = zeros(S,1);
        Num_initialization = (max(L,K)*10); % Số lần khởi tạo ngẫu nhiên ban đầu
        
        %% VÒNG LẶP MONTE CARLO
        for jj = 1:MonteCarlo
            % Tạo kênh vật lý G và kênh mục tiêu H_true
            G_independent = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
            G = sqrt(pathloss)*(Corr_R)^(1/2)*G_independent*(Corr_T)^(1/2);
            [~, G_svd, ~] = svd(G); 
            H_true = G_svd(1:S,1:S); 
            H_true_vec = H_true(:); % Fix lỗi vec
            Norm_H = norm(H_true_vec)^2;
            
            % Water-filling
            h_diag = diag(H_true);
            if S == 1
                PA_WF = Pt;
            else
                if exist('WF','file'), [PA_WF] = WF(Pt, Sigma2, h_diag); else, PA_WF = (Pt/S)*ones(S,1); end
            end
            
            % --- PHẦN 1: TÌM ĐIỂM KHỞI TẠO TỐT NHẤT (Initialization) ---
            Error_old_set = zeros(Num_initialization,1);
            phase_transmit_set = zeros(M,L,Num_initialization);
            phase_receive_set = zeros(N,K,Num_initialization);
            
            for tt = 1:Num_initialization
                phase_transmit = randn(M,L)+1i*randn(M,L); phase_transmit = phase_transmit./abs(phase_transmit);
                phase_receive = randn(N,K)+1i*randn(N,K); phase_receive = phase_receive./abs(phase_receive);
                
                % Tính P
                P = diag(phase_transmit(:,1))*W_T_1;
                for l=1:L-1, P = diag(phase_transmit(:,l+1))*W_T*P; end
                % Tính Q
                Q = U_R_1*diag(phase_receive(:,1));
                for k=1:K-1, Q = Q*U_R*diag(phase_receive(:,k+1)); end
                
                H_SIM = Q*G*P; H_SIM_vec = H_SIM(:);
                Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec;
                Error_old_set(tt) = norm(Factor*H_SIM_vec-H_true_vec)^2/Norm_H;
                phase_transmit_set(:,:,tt) = phase_transmit;
                phase_receive_set(:,:,tt) = phase_receive;
            end
            
            % Chọn bộ pha tốt nhất
            [~,d_max] = min(Error_old_set);
            Error_old = Error_old_set(d_max);
            phase_transmit = phase_transmit_set(:,:,d_max);
            phase_phase_transmit = angle(phase_transmit);
            phase_receive = phase_receive_set(:,:,d_max);
            phase_phase_receive = angle(phase_receive);
            
            % Tính lại H_SIM từ bộ pha tốt nhất
            P = diag(phase_transmit(:,1))*W_T_1;
            for l=1:L-1, P = diag(phase_transmit(:,l+1))*W_T*P; end
            Q = U_R_1*diag(phase_receive(:,1));
            for k=1:K-1, Q = Q*U_R*diag(phase_receive(:,k+1)); end
            H_SIM = Q*G*P; H_SIM_vec = H_SIM(:);
            Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec;

            % --- PHẦN 2: THUẬT TOÁN GRADIENT DESCENT (Iterative Algorithm) ---
            step = 0.1; 
            Error_new = 10000;
            
            while abs(Error_new-Error_old) >= Error_old * 0.001
                % 1. Tính đạo hàm TX
                for ll = 1:L
                    for mm = 1:M
                        X_left = W_T_1;
                        for ll_left = 1:ll-1, X_left = W_T*diag(phase_transmit(:,ll_left))*X_left; end
                        X_right = Q*G;
                        for ll_right = 1:(L-ll), X_right = X_right*diag(phase_transmit(:,L+1-ll_right))*W_T; end
                        for ss1 = 1:S
                            temp1 = X_right(:,mm)*X_left(mm,ss1);
                            Temp1(ss1) = 2*imag((Factor*phase_transmit(mm,ll)*temp1)'*(Factor*H_SIM(:,ss1)-H_true(:,ss1)));
                        end
                        Derivative_transmit_phase_shift(mm,ll) = sum(Temp1);
                    end
                end
                
                % 2. Tính đạo hàm RX
                for kk = 1:K
                    for nn = 1:N
                        Y_left = U_R_1;
                        for kk_left = 1:kk-1, Y_left = Y_left*diag(phase_receive(:,kk_left))*U_R; end
                        Y_right = G*P;
                        for kk_right = 1:(K-kk), Y_right = U_R*diag(phase_receive(:,K+1-kk_right))*Y_right; end
                        for ss1 = 1:S
                            Y = Y_left(ss1,nn)*Y_right(nn,:);
                            Temp2(ss1) = 2*imag((Factor*H_SIM(ss1,:)-H_true(ss1,:))*(Factor*phase_receive(nn,kk)*Y)');
                        end
                        Derivative_receive_phase_shift(nn,kk) = sum(Temp2);
                    end
                end
                
                % 3. Cập nhật pha
                Derivative_transmit_phase_shift = pi*Derivative_transmit_phase_shift/max(max(Derivative_transmit_phase_shift));
                phase_phase_transmit = phase_phase_transmit-step*Derivative_transmit_phase_shift;
                phase_transmit = exp(1i*phase_phase_transmit);
                
                Derivative_receive_phase_shift = pi*Derivative_receive_phase_shift/max(max(Derivative_receive_phase_shift));
                phase_phase_receive = phase_phase_receive-step*Derivative_receive_phase_shift;
                phase_receive = exp(1i*phase_phase_receive);
                
                step = step*0.5; % Giảm step size
                
                % 4. Tính lại Error
                P = diag(phase_transmit(:,1))*W_T_1;
                for l=1:L-1, P = diag(phase_transmit(:,l+1))*W_T*P; end
                Q = U_R_1*diag(phase_receive(:,1));
                for k=1:K-1, Q = Q*U_R*diag(phase_receive(:,k+1)); end
                H_SIM = Q*G*P; H_SIM_vec = H_SIM(:);
                Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec;
                
                Error_old = Error_new;
                Error_new = norm(Factor*H_SIM-H_true)^2/Norm_H;
            end
            
            % --- PHẦN 3: LƯU KẾT QUẢ ---
            NMSE(jj) = Error_new;
            C_single_stream = zeros(S,1);
            for pp = 1:S
                C_single_stream(pp) = log2(1+PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2/ ...
                    (Sigma2+(abs(Factor*H_SIM(pp,:)).^2*PA_WF-PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2)));
            end
            Capacity(jj) = sum(C_single_stream);
        end
        
        Capacity_average(ii) = mean(Capacity);
        NMSE_average(ii) = mean(NMSE);
        fprintf(' Done. Cap=%.2f\n', Capacity_average(ii));
        toc
    end
    
    % LƯU KẾT QUẢ RIÊNG CHO TỪNG K
    FileName_Cap = sprintf('Results/Benchmark_Capacity_K%d.mat', K);
    FileName_NMSE = sprintf('Results/Benchmark_NMSE_K%d.mat', K);
    
    if ~exist('Results', 'dir'), mkdir('Results'); end
    save(FileName_Cap, 'Capacity_average');
    save(FileName_NMSE, 'NMSE_average');
end
fprintf('\n--- DA CHAY XONG TOAN BO BENCHMARK (K=1,2,5,10) ---\n');