function [Best_Phi_Tx, Best_Phi_Rx, Best_Cost, Convergence_Curve] = PSO_Solver_SIM(SysModel, Init_Phi_Tx, Init_Phi_Rx)
% PSO_SOLVER_SIM: Phiên bản hỗ trợ Warm Start (Nhận 3 tham số đầu vào)

    %% 1. THIẾT LẬP THAM SỐ PSO
    n_particles = 40;       % Số hạt
    max_iter = 100;         % Số vòng lặp
    
    % Tham số quán tính và học tập
    w = 0.9;                
    w_damp = 0.98;          
    c1 = 1.49;              
    c2 = 1.49;              
    
    var_min = 0;
    var_max = 2*pi;
    
    %% 2. CHUẨN BỊ DỮ LIỆU
    M = SysModel.M; L = SysModel.L;
    N = SysModel.N; K = SysModel.K;
    
    dim_Tx = M * L;
    dim_Rx = N * K;
    n_vars = dim_Tx + dim_Rx;

    %% 3. KHỞI TẠO ĐÀN HẠT
    Particle_Pos = var_min + (var_max - var_min) * rand(n_particles, n_vars);
    Particle_Vel = zeros(n_particles, n_vars);
    Particle_Cost = inf(n_particles, 1);
    
    % --- WARM START: CẤY GHÉP CON ĐẦU ĐÀN (QUAN TRỌNG) ---
    % Kiểm tra xem có dữ liệu Warm Start được truyền vào không
    if nargin > 1 && ~isempty(Init_Phi_Tx) && ~isempty(Init_Phi_Rx)
        % Chuyển ma trận pha phức thành góc (Angle)
        angle_Tx = angle(Init_Phi_Tx);
        angle_Tx(angle_Tx < 0) = angle_Tx(angle_Tx < 0) + 2*pi;
        
        angle_Rx = angle(Init_Phi_Rx);
        angle_Rx(angle_Rx < 0) = angle_Rx(angle_Rx < 0) + 2*pi;
        
        % Ghép thành vector vị trí cho hạt số 1
        Warm_Start_Pos = [angle_Tx(:); angle_Rx(:)]';
        
        % Gán cho hạt đầu tiên để dẫn đường
        Particle_Pos(1, :) = Warm_Start_Pos;
    end
    % -----------------------------------------------------
    
    pBest_Pos = Particle_Pos;
    pBest_Cost = Particle_Cost;
    
    gBest_Pos = zeros(1, n_vars);
    gBest_Cost = inf;
    
    Convergence_Curve = zeros(max_iter, 1);

    %% 4. VÒNG LẶP PSO
    for it = 1:max_iter
        
        for i = 1:n_particles
            % Decode vị trí thành ma trận pha
            current_vec = Particle_Pos(i, :);
            
            vec_Tx = current_vec(1:dim_Tx);
            Phi_Tx = exp(1i * reshape(vec_Tx, [M, L]));
            
            vec_Rx = current_vec(dim_Tx+1:end);
            Phi_Rx = exp(1i * reshape(vec_Rx, [N, K]));
            
            % Tính P
            P = diag(Phi_Tx(:,1)) * SysModel.W_T_1;
            for l = 1:L-1
                P = diag(Phi_Tx(:,l+1)) * SysModel.W_T * P;
            end
            
            % Tính Q
            Q = SysModel.U_R_1 * diag(Phi_Rx(:,1));
            for k = 1:K-1
                Q = Q * SysModel.U_R * diag(Phi_Rx(:,k+1));
            end
            
            H_SIM = Q * SysModel.G * P;
            
            % Tính MSE
            H_SIM_vec = H_SIM(:);
            H_true_vec = SysModel.H_true(:);
            Factor = (H_SIM_vec' * H_SIM_vec) \ (H_SIM_vec' * H_true_vec);
            Current_MSE = norm(Factor * H_SIM_vec - H_true_vec)^2 / SysModel.Norm_H;
            
            % Cập nhật pBest
            if Current_MSE < pBest_Cost(i)
                pBest_Cost(i) = Current_MSE;
                pBest_Pos(i, :) = Particle_Pos(i, :);
                
                % Cập nhật gBest
                if pBest_Cost(i) < gBest_Cost
                    gBest_Cost = pBest_Cost(i);
                    gBest_Pos = pBest_Pos(i, :);
                end
            end
        end
        
        % Cập nhật vận tốc và vị trí
        r1 = rand(n_particles, n_vars);
        r2 = rand(n_particles, n_vars);
        
        Particle_Vel = w * Particle_Vel ...
                     + c1 * r1 .* (pBest_Pos - Particle_Pos) ...
                     + c2 * r2 .* (repmat(gBest_Pos, n_particles, 1) - Particle_Pos);
        
        % Giới hạn vận tốc
        max_vel = 0.2 * (var_max - var_min);
        Particle_Vel = max(min(Particle_Vel, max_vel), -max_vel);
                 
        Particle_Pos = Particle_Pos + Particle_Vel;
        
        % Xử lý biên
        Particle_Pos = max(min(Particle_Pos, var_max), var_min);
        
        Convergence_Curve(it) = gBest_Cost;
        w = w * w_damp; % Giảm quán tính
    end
    
    %% 5. TRẢ VỀ KẾT QUẢ
    final_vec = gBest_Pos;
    vec_Tx = final_vec(1:dim_Tx);
    Best_Phi_Tx = exp(1i * reshape(vec_Tx, [M, L]));
    vec_Rx = final_vec(dim_Tx+1:end);
    Best_Phi_Rx = exp(1i * reshape(vec_Rx, [N, K]));
    Best_Cost = gBest_Cost;
end