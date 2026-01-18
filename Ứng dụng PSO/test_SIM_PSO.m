% test_SIM_PSO_All_K.m
% Code này chạy thuật toán PSO (có Warm Start) với K = 1, 2, 5, 10
clc; clearvars; close all;

%% 1. CẤU HÌNH
K_values = [1, 2, 5, 10];
MonteCarlo = 20; % Số lần lặp (Tăng lên 20-50 để đường vẽ mượt hơn)
Max_L = 10;      % Số lớp phát

% Thông số vật lý
Thickness = 0.05; Pt = 10^(20/10); Sigma2 = 10^(-110/10);
c = 3*10^8; f0 = 28*10^9; lambda = c/f0;
N_max = 10; PL = -20*log10(4*pi/lambda)-35*log10(250); pathloss = 10^(PL/10);
M = 100; N = 100; d_element_spacing = lambda/2; S = 4;

%% 2. VÒNG LẶP CHÍNH
for k_idx = 1:length(K_values)
    K = K_values(k_idx);
    fprintf('\n>>> DANG CHAY PSO (Warm Start) VOI K = %d <<<\n', K);
    
    Capacity_average = zeros(Max_L,1);
    NMSE_average = zeros(Max_L,1);
    
    for ii = 1:Max_L
        L = ii;
        fprintf('   - Layer L=%d...', L);
        tic
        
        % --- TẠO KÊNH TRUYỀN (GIỐNG CODE GỐC) ---
        W_T = zeros(M,M); Corr_T = zeros(M,M); U_R = zeros(N,N); Corr_R = zeros(N,N);
        W_T_1 = zeros(M,S); U_R_1 = zeros(S,N);
        
        % Tính toán các ma trận kênh (Đã rút gọn)
        for mm1=1:M, for mm2=1:M, d_temp=sqrt((mod(mm1-1,10)-mod(mm2-1,10))^2+(ceil(mm1/10)-ceil(mm2/10))^2)*lambda/2; W_T(mm2,mm1)=lambda/4/pi/sqrt((Thickness/L)^2+d_temp^2)*exp(-1i*2*pi*sqrt((Thickness/L)^2+d_temp^2)/lambda); if d_temp==0, Corr_T(mm2,mm1)=1; else, Corr_T(mm2,mm1)=sin(pi*2*d_temp/lambda)/(pi*2*d_temp/lambda); end; end; end
        for nn1=1:N, for nn2=1:N, d_temp=sqrt((mod(nn1-1,10)-mod(nn2-1,10))^2+(ceil(nn1/10)-ceil(nn2/10))^2)*lambda/2; U_R(nn2,nn1)=lambda/4/pi/sqrt((Thickness/K)^2+d_temp^2)*exp(-1i*2*pi*sqrt((Thickness/K)^2+d_temp^2)/lambda); if d_temp==0, Corr_R(nn2,nn1)=1; else, Corr_R(nn2,nn1)=sin(pi*2*d_temp/lambda)/(pi*2*d_temp/lambda); end; end; end
        for mm=1:M, for nn=1:S, d_tr=sqrt((Thickness/L)^2+((mod(mm-1,10)+1-(1+10)/2)*lambda/2)^2+((ceil(mm/10)-(1+10)/2)*lambda/2-(nn-(1+S)/2)*lambda/2)^2); W_T_1(mm,nn)=lambda/4/pi/d_tr*exp(-1i*2*pi*d_tr/lambda); end; end
        for mm=1:N, for nn=1:S, d_re=sqrt((Thickness/K)^2+((mod(mm-1,10)+1-(1+10)/2)*lambda/2)^2+((ceil(mm/10)-(1+10)/2)*lambda/2-(nn-(1+S)/2)*lambda/2)^2); U_R_1(nn,mm)=lambda/4/pi/d_re*exp(-1i*2*pi*d_re/lambda); end; end

        Capacity = zeros(MonteCarlo,1); 
        NMSE = zeros(MonteCarlo,1);
        
        for jj = 1:MonteCarlo
            % Tạo kênh vật lý
            G_independent = sqrt(1/2)*(randn(N,M)+1i*randn(N,M));
            G = sqrt(pathloss)*(Corr_R)^(1/2)*G_independent*(Corr_T)^(1/2);
            [~, G_svd, ~] = svd(G); H_true = G_svd(1:S,1:S); H_true_vec = H_true(:); Norm_H = norm(H_true_vec)^2;
            
            % Water-filling đơn giản
            PA_WF = (Pt/S)*ones(S,1);

            % --- GIAI ĐOẠN 1: RANDOM SEARCH (TÌM ĐIỂM KHỞI ĐẦU) ---
            % Đây là phần còn thiếu ở code trước. Nó giúp PSO không bị "ngáo".
            Num_Init = 20; % Số lần thử ngẫu nhiên (chạy nhanh)
            Best_Rand_Error = inf;
            Best_Rand_Tx = []; Best_Rand_Rx = [];
            
            for tt = 1:Num_Init
                % Random pha
                pt = randn(M,L)+1i*randn(M,L); pt = pt./abs(pt);
                pr = randn(N,K)+1i*randn(N,K); pr = pr./abs(pr);
                
                % Tính nhanh MSE
                P_t = diag(pt(:,1))*W_T_1; for l=1:L-1, P_t = diag(pt(:,l+1))*W_T*P_t; end
                Q_t = U_R_1*diag(pr(:,1)); for k=1:K-1, Q_t = Q_t*U_R*diag(pr(:,k+1)); end
                H_S = Q_t*G*P_t; H_S_v = H_S(:);
                Fact = (H_S_v'*H_S_v)\H_S_v'*H_true_vec;
                Err = norm(Fact*H_S_v-H_true_vec)^2/Norm_H;
                
                if Err < Best_Rand_Error
                    Best_Rand_Error = Err;
                    Best_Rand_Tx = pt; Best_Rand_Rx = pr;
                end
            end
            
            % --- GIAI ĐOẠN 2: PSO (TỪ ĐIỂM TỐT NHẤT CỦA GIAI ĐOẠN 1) ---
            % Đóng gói struct
            SysModel.M=M; SysModel.N=N; SysModel.L=L; SysModel.K=K;
            SysModel.W_T=W_T; SysModel.W_T_1=W_T_1; SysModel.U_R=U_R; SysModel.U_R_1=U_R_1;
            SysModel.G=G; SysModel.H_true=H_true; SysModel.Norm_H=Norm_H;
            
            % Gọi PSO với Warm Start inputs
            [phase_transmit, phase_receive, Error_final, ~] = PSO_Solver_SIM(SysModel, Best_Rand_Tx, Best_Rand_Rx);
            
            % Lưu kết quả
            NMSE(jj) = Error_final;
            
            % Tính dung lượng cuối cùng
            P = diag(phase_transmit(:,1))*W_T_1; for l=1:L-1, P = diag(phase_transmit(:,l+1))*W_T*P; end
            Q = U_R_1*diag(phase_receive(:,1)); for k=1:K-1, Q = Q*U_R*diag(phase_receive(:,k+1)); end
            H_SIM = Q*G*P; H_SIM_vec = H_SIM(:); Factor = (H_SIM_vec'*H_SIM_vec)\H_SIM_vec'*H_true_vec; 
            
            C_stream = zeros(S,1);
            for pp=1:S, C_stream(pp)=log2(1+PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2/(Sigma2+(abs(Factor*H_SIM(pp,:)).^2*PA_WF-PA_WF(pp)*abs(Factor*H_SIM(pp,pp))^2))); end
            Capacity(jj) = sum(C_stream);
        end
        
        Capacity_average(ii) = mean(Capacity); 
        NMSE_average(ii) = mean(NMSE);
        fprintf(' Done. Cap=%.2f\n', Capacity_average(ii));
    end
    
    % LƯU FILE KẾT QUẢ
    if ~exist('Results', 'dir'), mkdir('Results'); end
    save(sprintf('Results/PSO_Capacity_K%d.mat', K), 'Capacity_average');
    save(sprintf('Results/PSO_NMSE_K%d.mat', K), 'NMSE_average');
end
fprintf('\n>>> DA HOAN THANH TAT CA! <<<\n');