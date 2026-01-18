#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <limits>
#include <random>
#include <chrono>
#include <iomanip>

using namespace std;

// Hằng số cần thiết
const double PI = acos(-1.0);
const double E = exp(1.0); // Hằng số e cho hàm Ackley

// --- Cấu hình Tham số Toàn cục (Global PSO Parameters) ---
// 1. Tham số Không gian Tìm kiếm
const double X_MIN = -5.0;  // Giới hạn dưới của không gian tìm kiếm
const double X_MAX = 5.0;   // Giới hạn trên của không gian tìm kiếm
const double VEL_MAX = (X_MAX - X_MIN); // Vận tốc tối đa ban đầu

// 2. Tham số Bầy đàn
const int POP_SIZE = 10;    // Kích thước bầy đàn (10 hạt)
const int MAX_ITER = 200;   // Số vòng lặp tối đa (200)

// 3. Tham số PSO
const double W = 0.7;       // Trọng số quán tính (Inertia Weight)
const double C1 = 1.5;      // Hệ số kinh nghiệm cá nhân (Cognitive factor)
const double C2 = 1.5;      // Hệ số kinh nghiệm tập thể (Social factor)

// Khởi tạo bộ sinh số ngẫu nhiên cho toàn bộ chương trình
mt19937_64 rng(chrono::high_resolution_clock::now().time_since_epoch().count());
// Phân phối chuẩn hóa (0.0 đến 1.0)
uniform_real_distribution<double> unif01(0.0, 1.0); 

// ======================================================================
// HÀM TIỆN ÍCH: Sinh số ngẫu nhiên trong khoảng [min, max]
// ======================================================================
double random_range(double min_val, double max_val) {
    // Công thức ánh xạ tuyến tính đơn giản: min + (số ngẫu nhiên [0, 1]) * (khoảng cách)
    return min_val + unif01(rng) * (max_val - min_val);
}

// Định nghĩa cấu trúc Hạt (Particle) trong bầy đàn
struct Particle {
    array<double, 2> pos;       // [x, y] - Vị trí hiện tại của hạt
    array<double, 2> vel;       // [vx, vy] - Vận tốc hiện tại của hạt
    array<double, 2> pbest;     // Vị trí tốt nhất cá nhân từng đạt được (Pbest)
    double pbest_fit;           // Giá trị fitness (hàm mục tiêu) tại Pbest
};

// Hàm mục tiêu: Ackley Function (tìm giá trị tối thiểu: f(0, 0) = 0)
double ackley(double x, double y) {
    return -20.0 * exp(-0.2 * sqrt(0.5 * (x*x + y*y)))
           - exp(0.5 * (cos(2.0 * PI * x) + cos(2.0 * PI * y)))
           + E + 20.0;
}

// ----------------------------------------------------------------------
// Hàm chính thực thi thuật toán PSO
// ----------------------------------------------------------------------
int main() {
    
    // Khai báo bầy đàn và biến lưu trữ kết quả toàn cục
    vector<Particle> swarm(POP_SIZE);
    array<double, 2> gbest_pos = {0.0, 0.0};
    double gbest_fit = numeric_limits<double>::infinity(); // Khởi tạo bằng giá trị lớn nhất

    // ------------------------------------------------------------------
    // 1. Khởi tạo Bầy đàn (Initialize Swarm)
    // ------------------------------------------------------------------
    cout << fixed << setprecision(6);
    cout << "=== KHOI TAO VI TRI NGẪU NHIÊN CHO " << POP_SIZE << " HAT ===\n";

    for (int i = 0; i < POP_SIZE; i++) {
        Particle &p = swarm[i];
        // sinh vị trí và vận tốc ngẫu nhiên
        p.pos = {random_range(X_MIN, X_MAX), random_range(X_MIN, X_MAX)};
        p.vel = {random_range(-VEL_MAX, VEL_MAX), random_range(-VEL_MAX, VEL_MAX)};
        
        // Thiết lập Pbest ban đầu
        p.pbest = p.pos;
        p.pbest_fit = ackley(p.pos[0], p.pos[1]);

        // Cập nhật Gbest ban đầu
        if (p.pbest_fit < gbest_fit) {
            gbest_fit = p.pbest_fit;
            gbest_pos = p.pbest;
        }

        // In ra vị trí khởi tạo
        cout << "Hat " << setw(2) << i + 1 << ": ";
        cout << "(" << p.pos[0] << ", " << p.pos[1] << ")\n";
    }

    cout << "\n=======================================================\n";

    // ------------------------------------------------------------------
    // 2. Vòng lặp chính PSO (Main Optimization Loop)
    // ------------------------------------------------------------------
    for (int iter = 1; iter <= MAX_ITER; iter++) {
        // Cập nhật từng hạt
        for (auto &p : swarm) {
            
            // Khai báo các biến ngẫu nhiên CỤC BỘ (r1, r2) đơn giản
            const double r1 = unif01(rng); 
            const double r2 = unif01(rng);

            // Cập nhật vận tốc và vị trí cho từng chiều (d=0 là x, d=1 là y)
            for (int d = 0; d < 2; d++) {
                
                // Công thức Cập nhật Vận tốc (PSO Formula)
                p.vel[d] = W * p.vel[d] 
                         + C1 * r1 * (p.pbest[d] - p.pos[d]) 
                         + C2 * r2 * (gbest_pos[d] - p.pos[d]);

                // Cập nhật Vị trí
                p.pos[d] += p.vel[d];

                // Giới hạn vị trí trong không gian tìm kiếm [X_MIN, X_MAX]
                if (p.pos[d] < X_MIN) p.pos[d] = X_MIN;
                if (p.pos[d] > X_MAX) p.pos[d] = X_MAX;
            }

            // Đánh giá và Cập nhật Pbest
            const double current_fit = ackley(p.pos[0], p.pos[1]);
            
            if (current_fit < p.pbest_fit) {
                p.pbest_fit = current_fit;
                p.pbest = p.pos;

                // Cập nhật Gbest (Tốt nhất toàn cục)
                if (current_fit < gbest_fit) {
                    gbest_fit = current_fit;
                    gbest_pos = p.pos;
                }
            }
        }
        

        if (iter % 10 == 0 || iter == 1) {
            cout << "Iter " << setw(4) << iter
                 << " | Best fitness: " << setw(12) << gbest_fit
                 << " | Best pos: (" << gbest_pos[0] << ", " << gbest_pos[1] << ")\n";
        }
    }

    // ------------------------------------------------------------------
    // 3. Kết quả Cuối cùng
    // ------------------------------------------------------------------
    cout << "\n=== KET QUA CUOI CUNG SAU " << MAX_ITER << " LAN LAP ===\n";
    cout << "Best fitness (Gia tri toi uu) = " << setprecision(8) << gbest_fit << "\n";
    cout << "Position (Vi tri): (" << gbest_pos[0] << ", " << gbest_pos[1] << ")\n";

    return 0;
}