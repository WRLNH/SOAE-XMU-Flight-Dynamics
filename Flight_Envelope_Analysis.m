m = 17; cbar = 0.423; S = 1.3536; b = 3.2; g = 9.81;
W = m * g;
IDX_alpha_deg = [-4; -2; 0; 2; 4; 8; 12; 16; 20];
TBL_C_L = [-0.219; -0.04; 0.139; 0.299; 0.455; 0.766; 1.083; 1.409; 1.743];
TBL_C_D = [0.026; 0.024; 0.024; 0.028; 0.036; 0.061; 0.102; 0.141; 0.173];
alpha_deg = -4;
V = zeros(1, 200);
H = zeros(1, 6);
Ma = zeros(1, 200);
D = zeros(6, 200);
for j = 1000 : 1000 : 6000
    for i = 1 : 200
        V(i) = i;
        H(j / 1000) = j;
        [rho, Ma(i)] = UAV_density(H(j / 1000), V(i));
        C_L = 2 * W / (rho * V(i)^2 * S);
        alpha_deg = -4;
        for idx = -4 : 0.0001 : 20
            C_L_0 = interp1d(TBL_C_L, IDX_alpha_deg, alpha_deg);
            if C_L_0 >= C_L
                break;
            end
            alpha_deg = alpha_deg + 0.0001;
        end
        C_D = interp1d(TBL_C_D, IDX_alpha_deg, alpha_deg);
        D(j / 1000, i) = C_D * rho * V(i)^2 * S / 2;
    end
end
figure;
plot(Ma, D);
legend('H=1000', 'H=2000', 'H=3000', 'H=4000', 'H=5000', 'H=6000');
grid on;
