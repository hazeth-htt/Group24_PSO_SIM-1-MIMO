function [Best_Phi_Tx, Best_Phi_Rx, Best_Cost, Convergence_Curve] = PSO_Solver_DEPSO(SysModel, Init_Phi_Tx, Init_Phi_Rx)
% PSO_SOLVER_DEPSO: Biến thể Differential Evolution + PSO
% Hàm này độc lập với PSO thường để bạn có thể so sánh cả hai.

    %% 1. THIẾT LẬP THAM SỐ DEPSO
    n_particles = 40;       
    max_iter = 100;         
    
    % Tham số PSO
    w = 0.73;               
    c1 = 1.5; c2 = 1.5;     
    
    % Tham số DE
    beta_min = 0.2; beta_max = 0.8; % Hệ số khuếch đại
    CR = 0.5;                       % Xác suất lai ghép (Crossover)
    
    var_min = 0; var_max = 2*pi;
    max_vel = 0.2 * (var_max - var_min);
    
    %% 2. CHUẨN BỊ
    M = SysModel.M; L = SysModel.L;
    N = SysModel.N; K = SysModel.K;
    dim_Tx = M * L; dim_Rx = N * K; n_vars = dim_Tx + dim_Rx;

    %% 3. KHỞI TẠO & WARM START
    Particle_Pos = var_min + (var_max - var_min) * rand(n_particles, n_vars);
    Particle_Vel = zeros(n_particles, n_vars);
    Particle_Cost = inf(n_particles, 1);
    
    if nargin > 1 && ~isempty(Init_Phi_Tx) && ~isempty(Init_Phi_Rx)
        angle_Tx = angle(Init_Phi_Tx); angle_Tx(angle_Tx < 0) = angle_Tx(angle_Tx < 0) + 2*pi;
        angle_Rx = angle(Init_Phi_Rx); angle_Rx(angle_Rx < 0) = angle_Rx(angle_Rx < 0) + 2*pi;
        Warm_Start_Pos = [angle_Tx(:); angle_Rx(:)]';
        Particle_Pos(1, :) = Warm_Start_Pos;
    end
    
    pBest_Pos = Particle_Pos; pBest_Cost = Particle_Cost;
    gBest_Pos = zeros(1, n_vars); gBest_Cost = inf;
    Convergence_Curve = zeros(max_iter, 1);

    % Đánh giá ban đầu
    for i = 1:n_particles
        Cost = CostFunc(Particle_Pos(i,:), SysModel, M, L, N, K);
        Particle_Cost(i) = Cost;
        if Cost < pBest_Cost(i)
            pBest_Cost(i) = Cost; pBest_Pos(i,:) = Particle_Pos(i,:);
            if pBest_Cost(i) < gBest_Cost, gBest_Cost = pBest_Cost(i); gBest_Pos = pBest_Pos(i,:); end
        end
    end

    %% 4. VÒNG LẶP CHÍNH
    for it = 1:max_iter
        beta = beta_min + (beta_max - beta_min) * rand(); % Adaptive Beta
        
        for i = 1:n_particles
            % --- 1. Cập nhật PSO ---
            r1 = rand(1, n_vars); r2 = rand(1, n_vars);
            Particle_Vel(i,:) = w*Particle_Vel(i,:) + c1*r1.*(pBest_Pos(i,:) - Particle_Pos(i,:)) + c2*r2.*(gBest_Pos - Particle_Pos(i,:));
            Particle_Vel(i,:) = max(min(Particle_Vel(i,:), max_vel), -max_vel);
            New_Pos_PSO = Particle_Pos(i,:) + Particle_Vel(i,:);
            New_Pos_PSO = max(min(New_Pos_PSO, var_max), var_min);
            
            % --- 2. Lai ghép DE ---
            idxs = randperm(n_particles, 3);
            while any(idxs == i), idxs = randperm(n_particles, 3); end
            r1=idxs(1); r2=idxs(2); r3=idxs(3);
            
            Mutant = pBest_Pos(r1,:) + beta * (pBest_Pos(r2,:) - pBest_Pos(r3,:));
            Mutant = max(min(Mutant, var_max), var_min);
            
            Trial_Pos = New_Pos_PSO;
            j_rand = randi(n_vars);
            mask = (rand(1, n_vars) <= CR); mask(j_rand) = 1;
            Trial_Pos(mask) = Mutant(mask);
            
            % --- 3. Chọn lọc ---
            Cost_Trial = CostFunc(Trial_Pos, SysModel, M, L, N, K);
            if Cost_Trial < Particle_Cost(i)
                Particle_Pos(i,:) = Trial_Pos;
                Particle_Cost(i) = Cost_Trial;
                if Particle_Cost(i) < pBest_Cost(i)
                    pBest_Cost(i) = Particle_Cost(i); pBest_Pos(i,:) = Particle_Pos(i,:);
                    if pBest_Cost(i) < gBest_Cost, gBest_Cost = pBest_Cost(i); gBest_Pos = pBest_Pos(i,:); end
                end
            else
                Particle_Pos(i,:) = New_Pos_PSO; % Giữ hướng PSO nếu DE không tốt hơn
            end
        end
        Convergence_Curve(it) = gBest_Cost;
        w = w * 0.98;
    end
    
    %% 5. KẾT QUẢ
    final_vec = gBest_Pos;
    dim_Tx = M*L;
    vec_Tx = final_vec(1:dim_Tx); Best_Phi_Tx = exp(1i * reshape(vec_Tx, [M, L]));
    vec_Rx = final_vec(dim_Tx+1:end); Best_Phi_Rx = exp(1i * reshape(vec_Rx, [N, K]));
    Best_Cost = gBest_Cost;
end

function C = CostFunc(v, S, M, L, N, K)
    dt=M*L; vt=v(1:dt); Pt=exp(1i*reshape(vt,[M,L]));
    vr=v(dt+1:end); Pr=exp(1i*reshape(vr,[N,K]));
    
    P=diag(Pt(:,1))*S.W_T_1; for l=1:L-1, P=diag(Pt(:,l+1))*S.W_T*P; end
    Q=S.U_R_1*diag(Pr(:,1)); for k=1:K-1, Q=Q*S.U_R*diag(Pr(:,k+1)); end
    
    H=Q*S.G*P; Hv=H(:); Ht=S.H_true(:);
    F=(Hv'*Hv)\Hv'*Ht;
    C=norm(F*Hv-Ht)^2/S.Norm_H;
end