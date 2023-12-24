V = 35; H = 1200; alpha_deg = 0.4495; theta_deg = alpha_deg; eng = 0.539353;
m = 17; Ix = 1.71; Iy = 5.54; Iz = 4.15; Izx = 0;
cbar = 0.423; SA = 1.3536; b = 3.2; g = 9.81;
Myrad2deg = 180 / pi;
T = eng / 100 * (m * g) / 4;
[ru, mach] = UAV_density(H, V);
q = (ru * V * V / 2);
alpha = alpha_deg / Myrad2deg; theta = theta_deg / Myrad2deg;

[CD0] = UAV_CD(alpha_deg);
[CD1] = UAV_CD(alpha_deg + 0.5);
[CD2] = UAV_CD(alpha_deg - 0.5);
CDalpha = (CD1 - CD2) * Myrad2deg;
[CL0, ~] = UAV_CL(alpha_deg);
[CL1, ~] = UAV_CL(alpha_deg + 0.5);
[CL2, CL_ele] = UAV_CL(alpha_deg - 0.5);
CLalpha = (CL1 - CL2) * Myrad2deg;

[CM0, ~, ~, ~] = UAV_CM(alpha_deg);
[CM1, ~, ~, ~] = UAV_CM(alpha_deg + 0.5);
[CM2, ~, CM_Q, CM_dalpha] = UAV_CM(alpha_deg - 0.5);
CMalpha = (CM1 - CM2) * Myrad2deg;

Xv = (-2 * CD0 * q * SA) / (m * V);
Xa = (-T * sin(alpha) - CDalpha * q * SA) / m;

Zv = (2 * CL0 * q * SA) / (m * V * V);
Za = (T * cos(alpha) + CLalpha * q * SA) / (m * V);

Mbar_V = 0;
Mbar_dalpha = CM_dalpha * (cbar / (2 * V)) * ((q * SA * cbar) / (Iy));
Mbar_alpha = CMalpha * ((q * SA * cbar) / (Iy));
Mbar_q = CM_Q * (cbar / (2 * V)) * ((q * SA * cbar) / (Iy));

Myalong = [Xv, Xa + g, 0, -g;
           -Zv, -Za, 1, 0;
           Mbar_V - Mbar_dalpha * Zv, Mbar_alpha - Mbar_dalpha * Za, Mbar_q + Mbar_dalpha, 0;
           0, 0, 1, 0];

Xdeltap = T * cos(alpha) / m;

Zdeltae = CL_ele * ((q * SA) / (m * V));
Zdeltap = T * sin(alpha) / (m * V);

CM_ele = -0.02052;
Mbar_deltae = CM_ele * (q * SA * cbar) / Iy;

Myblong = [0, Xdeltap;
           -Zdeltae, -Zdeltap;
           Mbar_deltae - Mbar_dalpha * Zdeltae, 0 - Mbar_dalpha * Zdeltap;
           0, 0];

Ybar = [0, 0, 0, 0, 0];
N = [0, 0, 0, 0, 0];
L = [0, 0, 0, 0, 0];
Nbar = [0, 0, 0, 0, 0];
Lbar = [0, 0, 0, 0, 0];
C_y = [-0.00909 * Myrad2deg, 0.03417, 0.42111, 0.000294, -0.002616];
C_n = [0.00235 * Myrad2deg, -0.01792, -0.15844, 0.000132, -0.00111];
C_l = [-0.00600 * Myrad2deg, -0.52568, 0.01832, -0.003618, -0.000144];

for i = 1:5
    Ybar(i) = (C_y(i) * q * SA) / (m * V);
    N(i) = C_n(i) * q * SA * b;
    L(i) = C_l(i) * q * SA * b;

    if i == 2 || i == 3
        Ybar(i) = Ybar(i) * b / (2 * V);
        N(i) = N(i) * b / (2 * V);
        L(i) = L(i) * b / (2 * V);
    end

end

for i = 1:5
    Lbar(i) = L(i) / Ix;
    Nbar(i) = N(i) / Iz;
end

Myalate = [Ybar(1), alpha + Ybar(2), Ybar(3) - 1, g * cos(theta) / V;
           Lbar(1), Lbar(2), Lbar(3), 0;
           Nbar(1), Nbar(2), Nbar(3), 0;
           0, 1, tan(theta), 0];

Myblate = [0, Ybar(5);
           Lbar(4), Lbar(5);
           Nbar(4), Nbar(5);
           0, 0];

CMele_over_CMalpha = CM_ele / CM0;

CNrud_over_CNbeta = -0.00111/0.00235;

omega_nsp = sqrt(- (Mbar_alpha + Mbar_q * Za));
xi_sp =- (Mbar_q + Mbar_dalpha - Za) / (2 * omega_nsp);
omega_np = sqrt(2) * g / V;
xi_p = 1 / (sqrt(2) * (CL0 / CD0));
omega_ndr = sqrt(Nbar(1) - Nbar(1) * Ybar(3) + Nbar(3) * Ybar(1));
xi_dr =- (Nbar(3) + Ybar(1)) / (2 * omega_ndr);
