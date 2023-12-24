%% 变量定义
m = 17; S = 1.3536; g = 9.81;
W = m * g;
i_alpha_deg = [-4; -2; 0; 2; 4; 8; 12; 16; 20];
TBL_C_L = [-0.219; -0.04; 0.139; 0.299; 0.455; 0.766; 1.083; 1.409; 1.743];
TBL_C_D = [0.026; 0.024; 0.024; 0.028; 0.036; 0.061; 0.102; 0.141; 0.173];
C_L_max = interp1d(TBL_C_L, i_alpha_deg, 4);

%% 平飞需求推力曲线
V = zeros(1, 200);
H = zeros(1, 6);
Ma_max = zeros(1, 200);
D = zeros(6, 200);

for j = 1:6

    for i = 1:200
        V(i) = i;
        H(j) = j * 1000;
        [rho, Ma_max(i)] = UAV_density(H(j), V(i));
        C_L = 2 * W / (rho * V(i) ^ 2 * S);

        for alpha_deg = -4:0.0001:20
            C_L_0 = interp1d(TBL_C_L, i_alpha_deg, alpha_deg);

            if C_L_0 >= C_L
                break;
            end

        end

        C_D = interp1d(TBL_C_D, i_alpha_deg, alpha_deg);
        D(j, i) = C_D * rho * V(i) ^ 2 * S / 2;
    end

end

figure;
plot(Ma_max, D);
xlabel("$Ma$", "Interpreter", "latex");
ylabel("$P$", "Interpreter", "latex");
legend('H=1000', 'H=2000', 'H=3000', 'H=4000', 'H=5000', 'H=6000');
title("平飞需求推力曲线");
grid on;

%% 最小速度曲线
rho = zeros(1, 6001);
H = zeros(1, 6001);
V_min = zeros(1, 6001);
Ma_min = zeros(1, 6001);

for i = 0:6000
    H(i + 1) = i;
    [rho(i + 1)] = My_density(H(i + 1));
    V_min(i + 1) = sqrt(2 * W / (C_L_max * rho(i + 1) * S));
    [~, Ma_min(i + 1)] = UAV_density(H(i + 1), V_min(i + 1));
end

xmin = max(Ma_min);
figure;
plot(Ma_min, H);
xlabel("$Ma$", "Interpreter", "latex");
ylabel("$H$", "Interpreter", "latex");
title("最小速度曲线");
grid on;

%% 最大速度曲线
C_D = interp1d(TBL_C_D, i_alpha_deg, -1.2849);
V_max = zeros(1, 6001);
Ma_max = zeros(1, 6001);

for i = 0:6000
    V_max(i + 1) = sqrt(2 * W / (C_D * rho(i + 1) * S));
    [~, Ma_max(i + 1)] = UAV_density(H(i + 1), V_max(i + 1));
end

xmax = max(Ma_max);
figure;
plot(Ma_max, H);
xlabel("$Ma$", "Interpreter", "latex");
ylabel("$H$", "Interpreter", "latex");
title("最大速度曲线");
grid on;

%% 飞行包线
figure;
plot([xmin, xmax], [6000, 6000], "b", Ma_min, H, "b", Ma_max, H, "b");
ylim([0, 6500]);
xlabel("$Ma$", "Interpreter", "latex");
ylabel("$H$", "Interpreter", "latex");
title("飞行包线");
grid on;
