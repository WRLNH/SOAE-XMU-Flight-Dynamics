V = 35; H = 1200; alpha_deg = 0.4495; eng = 0.5393;
% V = 25; H = 1000;alpha_deg = 2.7511; eng = 0.34986;
m = 17; Ix = 1.71; Iy = 5.54; Iz = 4.15;
cbar = 0.423; SA = 1.3536; b = 3.2; g = 9.81;
Myrad2deg = 180 / pi;
T = eng / 100 * (m * g) / 4;
[ru, mach] = UAV_density(H, V);
q = (ru * V * V / 2);
alpha = alpha_deg / Myrad2deg;

[CD0] = UAV_CD(alpha_deg);
[CD1] = UAV_CD(alpha_deg + 0.5);
[CD2] = UAV_CD(alpha_deg - 0.5);
CDalpha = (CD1 - CD2) * Myrad2deg;
[CL0, ~] = UAV_CL(alpha_deg);
[CL1, ~] = UAV_CL(alpha_deg + 0.5);
[CL2, CL_ele] = UAV_CL(alpha_deg - 0.5);
CLalpha = (CL1 - CL2) * Myrad2deg;
% CL_ele = CL_ele * Myrad2deg;
[CM0, ~, ~, ~] = UAV_CM(alpha_deg);
[CM1, ~, ~, ~] = UAV_CM(alpha_deg + 0.5);
[CM2, ~, CM_Q, CM_dalpha] = UAV_CM(alpha_deg - 0.5);
CMalpha = (CM1 - CM2) * Myrad2deg;
% CM_Q = CM_Q * Myrad2deg;
% CM_dalpha = CM_dalpha * Myrad2deg;

Xv = (-2 * CD0 * q * SA) / (m * V);
Xa = (-T * sin(alpha) - CDalpha * q * SA) / m;
along11 = Xv; along12 = Xa + g; along13 = 0; along14 = -g;

Zv = (2 * CL0 * q * SA) / (m * V * V);
Za = (T * cos(alpha) + CLalpha * q * SA) / (m * V);
along21 = -Zv; along22 = -Za; along23 = 1; along24 = 0;

Mbar_V = 0; % (2 * CM0 * q * SA * cbar) / (V * Iy);
Mbar_dalpha = CM_dalpha * (cbar / (2 * V)) * ((q * SA * cbar) / (Iy));
Mbar_alpha = CMalpha * ((q * SA * cbar) / (Iy));
Mbar_q = CM_Q * (cbar / (2 * V)) * ((q * SA * cbar) / (Iy));
along31 = Mbar_V - Mbar_dalpha * Zv; along32 = Mbar_alpha - Mbar_dalpha * Za; along33 = Mbar_q + Mbar_dalpha; along34 = 0;

along41 = 0; along42 = 0; along43 = 1; along44 = 0;

Myalong = [along11, along12, along13, along14;
           along21, along22, along23, along24;
           along31, along32, along33, along34;
           along41, along42, along43, along44];

Xdeltap = T * cos(alpha) / m;
blong11 = 0; blong12 = Xdeltap;

Zdeltae = CL_ele * ((q * SA) / (m * V));
Zdeltap = T * sin(alpha) / (m * V);
blong21 = -Zdeltae; blong22 = -Zdeltap;

CM_ele = -0.02052;
Mbar_deltae = CM_ele * (q * SA * cbar) / Iy;
blong31 = Mbar_deltae - Mbar_dalpha * Zdeltae; blong32 = 0 - Mbar_dalpha * Zdeltap;

blong41 = 0; blong42 = 0;

Myblong = [blong11, blong12;
           blong21, blong22;
           blong31, blong32;
           blong41, blong42];
