% test_SIM_DEPSO.m
% Chạy thuật toán biến thể DE-PSO (Phien ban Fix loi hoan chinh)
clc; clearvars; close all;

%% CẤU HÌNH
K_values = [1, 2, 5, 10];
MonteCarlo = 20; % So lan lap
Max_L = 10;
Thickness = 0.05; Pt = 10^(20/10); Sigma2 = 10^(-110/10);
c = 3*10^8; f0 = 28*10^9; lambda = c/f0;
N_max = 10; PL = -20*log10(4*pi/lambda)-35*log10(250); pathloss = 10^(PL/10);
M = 100; N = 100; d_element_spacing = lambda/2; S = 4;

for k_idx = 1:length(K_values)
    K = K_values(k_idx);
    fprintf('\n>>> DANG CHAY DEPSO (BIEN THE) VOI K = %d <<<\n', K);
    Capacity_average = zeros(Max_L,1); NMSE_average = zeros(Max_L,1);
    
    for ii = 1:Max_L
        L = ii; fprintf('   - Layer L=%d...', L); tic;
        
        % Tạo kênh (Rút gọn)
        W_T=zeros(M,M); W_T_1=zeros(M,S); U_R=zeros(N,N); U_R_1=zeros(S,N);
        d_layer_spacing_transmit = Thickness/L; d_layer_spacing_receive = Thickness/K;
        
        % Tinh toan ma tran kenh
        for mm1=1:M, for mm2=1:M, d_temp=sqrt((mod(mm1-1,10)-mod(mm2-1,10))^2+(ceil(mm1/10)-ceil(mm2/10))^2)*lambda/2; W_T(mm2,mm1)=lambda/4/pi/sqrt((Thickness/L)^2+d_temp^2)*exp(-1i*2*pi*sqrt((Thickness/L)^2+d_temp^2)/lambda); end; end
        for nn1=1:N, for nn2=1:N, d_temp=sqrt((mod(nn1-1,10)-mod(nn2-1,10))^2+(ceil(nn1/10)-ceil(nn2/10))^2)*lambda/2; U_R(nn2,nn1)=lambda/4/pi/sqrt((Thickness/K)^2+d_temp^2)*exp(-1i*2*pi*sqrt((Thickness/K)^2+d_temp^2)/lambda); end; end
        for mm=1:M, for nn=1:S, d_tr=sqrt((Thickness/L)^2+((mod(mm-1,10)+1-(1+10)/2)*lambda/2)^2+((ceil(mm/10)-(1+10)/2)*lambda/2-(nn-(1+S)/2)*lambda/2)^2); W_T_1(mm,nn)=lambda/4/pi/d_tr*exp(-1i*2*pi*d_tr/lambda); end; end
        for mm=1:N, for nn=1:S, d_re=sqrt((Thickness/K)^2+((mod(mm-1,10)+1-(1+10)/2)*lambda/2)^2+((ceil(mm/10)-(1+10)/2)*lambda/2-(nn-(1+S)/2)*lambda/2)^2); U_R_1(nn,mm)=lambda/4/pi/d_re*exp(-1i*2*pi*d_re/lambda); end; end
        
        Capacity = zeros(MonteCarlo,1); NMSE = zeros(MonteCarlo,1);
        
        for jj = 1:MonteCarlo
            G = sqrt(pathloss)*sqrt(1/2)*(randn(N,M)+1i*randn(N,M)); 
            [~, G_svd, ~] = svd(G); H_true = G_svd(1:S,1:S); Norm_H = norm(H_true(:))^2;
            
            % --- PHAN BO CONG SUAT (EQUAL POWER) ---
            % Su dung Equal Power de tranh loi Water-Filling
            PA_WF = (Pt/S)*ones(S,1);
            
            % --- WARM START RANDOM SEARCH ---
            Num_Init = 20; Best_Err = inf; BTx=[]; BRx=[];
            for tt=1:Num_Init
                pt=randn(M,L)+1i*randn(M,L); pt=pt./abs(pt);
                pr=randn(N,K)+1i*randn(N,K); pr=pr./abs(pr);
                Pt_mat=diag(pt(:,1))*W_T_1; for l=1:L-1, Pt_mat=diag(pt(:,l+1))*W_T*Pt_mat; end
                Qt_mat=U_R_1*diag(pr(:,1)); for k=1:K-1, Qt_mat=Qt_mat*U_R*diag(pr(:,k+1)); end
                Hs=Qt_mat*G*Pt_mat; Hv=Hs(:); Ht=H_true(:); F=(Hv'*Hv)\Hv'*Ht;
                if norm(F*Hv-Ht)^2/Norm_H < Best_Err, Best_Err=norm(F*Hv-Ht)^2/Norm_H; BTx=pt; BRx=pr; end
            end
            
            % --- GOI HAM DEPSO ---
            SM.M=M; SM.N=N; SM.L=L; SM.K=K; SM.W_T=W_T; SM.W_T_1=W_T_1; SM.U_R=U_R; SM.U_R_1=U_R_1; SM.G=G; SM.H_true=H_true; SM.Norm_H=Norm_H;
            [pT, pR, Err, ~] = PSO_Solver_DEPSO(SM, BTx, BRx);
            
            NMSE(jj) = Err;
            
            % Tinh toan lai kenh tong hop
            P=diag(pT(:,1))*W_T_1; for l=1:L-1, P=diag(pT(:,l+1))*W_T*P; end
            Q=U_R_1*diag(pR(:,1)); for k=1:K-1, Q=Q*U_R*diag(pR(:,k+1)); end
            H_SIM=Q*G*P; Hv=H_SIM(:); 
            
            % Tinh Factor (He so bu)
            F=(Hv'*Hv)\Hv'*H_true(:);
            
            % --- TINH DUNG LUONG (DA FIX LOI DIMENSION) ---
            C_s = zeros(S,1);
            PA_Vector = PA_WF(:); % Ep ve vector cot an toan
            
            for pp = 1:S
                % Xu ly Factor F (neu la vector hay scalar deu duoc)
                if length(F) > 1, F_val = F(pp); else, F_val = F; end
                
                H_row = H_SIM(pp, :); % Lay hang thu pp
                H_eff_row = F_val * H_row; % Kenh hieu dung (Vector Hang)
                
                % Tinh tong cong suat thu (Scalar)
                % (1xS) * (Sx1) -> Scalar
                Total_Power_Received = (abs(H_eff_row).^2) * PA_Vector;
                
                % Tinh tin hieu mong muon
                Signal_Power = PA_Vector(pp) * abs(H_eff_row(pp))^2;
                
                % Tinh Nhieu va Dung luong
                Interference_Plus_Noise = Sigma2 + (Total_Power_Received - Signal_Power);
                C_s(pp) = log2(1 + Signal_Power / Interference_Plus_Noise);
            end
            Capacity(jj) = sum(C_s);
            % ------------------------------------------------
        end
        Capacity_average(ii) = mean(Capacity); NMSE_average(ii) = mean(NMSE);
        toc
    end
    if ~exist('Results', 'dir'), mkdir('Results'); end
    save(sprintf('Results/DEPSO_Capacity_K%d.mat', K), 'Capacity_average');
    save(sprintf('Results/DEPSO_NMSE_K%d.mat', K), 'NMSE_average');
end