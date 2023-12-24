%% 变量定义
m = 17; S = 1.3536; g = 9.81; gamma = 14;
W = m * g;
IDX_alpha_deg = [-4; -2; 0; 2; 4; 8; 12; 16; 20];
TBL_C_L = [-0.219; -0.04; 0.139; 0.299; 0.455; 0.766; 1.083; 1.409; 1.743];
TBL_C_D = [0.026; 0.024; 0.024; 0.028; 0.036; 0.061; 0.102; 0.141; 0.173];

%% 不同高度的离地速度
L = W * cosd(gamma);
C_L_at_14 = interp1d(TBL_C_L, IDX_alpha_deg, 14);
H = zeros(1, 6001);
rho = zeros(1, 6001);
V_min = zeros(1, 6001);
Ma_min = zeros(1, 6001);

for i = 0:6000
    H(i + 1) = i;
    [rho(i + 1)] = My_density(H(i + 1));
    V_min(i + 1) = sqrt(2 * L / (C_L_at_14 * rho(i + 1) * S));
    [~, Ma_min(i + 1)] = UAV_density(H(i + 1), V_min(i + 1));
end

figure;
plot(Ma_min, H);
xlabel("$Ma$", "Interpreter", "latex");
ylabel("$H$", "Interpreter", "latex");
title("最小离地速度");
grid on;

%% 无动力下滑距离
H = 1000;
C_L = interp1d(TBL_C_L, IDX_alpha_deg, 4);
C_D = interp1d(TBL_C_D, IDX_alpha_deg, 4);
K = C_L / C_D;
R_d = H * K;
