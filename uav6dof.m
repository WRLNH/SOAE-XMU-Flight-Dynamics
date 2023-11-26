function [sys,x0,str,ts] = uav6dof(t,x,u,flag)
%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 12; %12个状态量
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 13; %自己定
sizes.NumInputs      = 4;  %4个输入量
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);
%
% initialize the initial conditions 提供初始状态量
%
rad2deg=57.295779513082320876798154814105;

psi_hmr=0.0;
if (psi_hmr> 180)   psi_hmr = psi_hmr-360.0; end
if (psi_hmr<-180)   psi_hmr = psi_hmr+360.0; end

Vt =30; alpha_deg=0.0; beta_deg =0.0; theta_deg=alpha_deg; H=1000;

x0  = [Vt; alpha_deg/rad2deg; beta_deg/rad2deg; 0; theta_deg/rad2deg;psi_hmr/rad2deg; 0;0;0; 0; 0; H];
%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
% ---状态量-------
Vt      = x(1);    phi   = x(4);    P = x(7);    PN = x(10);  
alpha   = x(2);    theta = x(5);    Q = x(8);    PE = x(11); 
beta    = x(3);    psi   = x(6);    R = x(9);     H = x(12);
rad2deg=57.295779513082320876798154814105;
alpha_deg = alpha * rad2deg;
% ---输入量-------
ele = u(1);      % elevator deflection angle (deg)
ail = u(2);      % aileron  deflection angle (deg)
rud = u(3);      % rudder   deflection angle (deg)
eng = u(4);   
% ---质量数据-------
mass = 17; Ix = 1.71; Iy = 5.54; Iz = 4.15; cbar = 0.423; SA = 1.3536; b = 3.2; g = 9.91;

% ---- 计算推力----
Pow =eng/100*(mass*g/3.0);
% ---- 计算空气密度和动压----
[ru,mach] = UAV_density(H,Vt);   % [air density] [mach number]
qs = SA*(ru*Vt*Vt/2);         % [Dynamic pressure](kg/m^2)
% ---- 计算气动力----
[CL0,CL_ele]=UAV_CL(alpha_deg);
L = qs*(CL0 + CL_ele*ele);
[CD0]=UAV_CD(alpha_deg);
D = qs*(CD0);
[CM0, CM_alpha] = UAV_CM(alpha_deg);

% ---- 基于气动数据计算气动受力，并以十二个方程求解状态量的微分量----

B = [cos(psi) * cos(theta), sin(psi) * cos(theta), -sin(theta);
     cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi), sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi), cos(theta) * sin(phi);
     cos(phi) * sin(theta) * cos(psi) + sin(psi) * sin(phi), sin(psi) * sin(theta) * sin(phi) - cos(psi) * sin(phi), cos(theta) * cos(phi)];
S = [cos(alpha) * cos(beta), sin(beta), sin(alpha) * cos(beta);
     -cos(alpha) * sin(beta), cos(beta), -sin(alpha) * sin(beta);
     -sin(alpha), 0, cos(alpha)];
ST = transpose(S);

CD = UAV_CD(alpha_deg);
CY = -0.00909 * beta;
CL = 0.00636 * ele + CL0 + CL_ele;
D = qs * CD;
Y = qs * CY;
L = qs * CL;

PQR = [P; Q; R];
U = Vt * cos(alpha) * cos(beta);
V = Vt * sin(beta);
W = Vt * sin(alpha) * cos(beta);
UVW = [U; V; W];
Fxyz = [Pow; 0; 0] + ST * [-D; Y; -L];
Fx = Fxyz(1); Fy = Fxyz(2); Fz = Fxyz(3);

dUVW = Fxyz / mass - cross(PQR, UVW) + g * [-sin(theta); sin(phi) * cos(theta); cos(phi) * cos(theta)];
dU = dUVW(1); dV = dUVW(2); dW = dUVW(3);
dVt = (U * dU + V * dV + W * dW) / Vt;
dbeta = (dV * Vt - V * dVt) / (Vt * Vt * cos(beta));
dalpha = (U * dW - W * dU) / (U * U + W * W);

dXYZ = B * UVW;

dX = dXYZ(1); dY = dXYZ(2); dZ = dXYZ(3);

CR = -0.006 * beta - 0.003618 * ail - 0.000144 * rud + b / (2 * Vt) * (-0.52568 * P + 0.01832 * R);
CM = CM0 + CM_alpha * alpha - 0.02052 * ele + cbar / (2 * Vt) * (-9.3136 * Q - 4.0258 * dalpha);
CN = 0.00235 * beta + 0.000132 * ail - 0.00111 * rud + b / (2 * Vt) * (-0.01792 * P - 0.15844 * R);
Lbar = qs * b * CR;
M = qs * cbar * CM;
N = qs * b * CN;

J = [Ix, 0, 0;
     0, Iy, 0;
     0, 0, Iz];

dPQR = J \ (-cross(PQR, (J * PQR)) + [Lbar; M; N]);
dP = dPQR(1); dQ = dPQR(2); dR = dPQR(3);

dPhi = P + (sin(theta) / cos(theta)) * (Q * sin(phi) + R * cos(phi));
dTheta = Q * cos(phi) - R * sin(phi);
dPsi = (Q * sin(phi) + R * cos(phi)) / cos(theta);

dPN = dX;
dPE = dY;
dH = -dZ;

% ---状态量的微分量-------
sys = [dVt;dalpha;dbeta;dPhi;dTheta;dPsi;dP;dQ;dR;dPN;dPE;dH;];
% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
% ---状态量-------
Vt     = x(1);    phi   = x(4);    P = x(7);    PN = x(10);  
alpha  = x(2);    theta = x(5);    Q = x(8);    PE = x(11);  
beta   = x(3);    psi   = x(6);    R = x(9);    H = x(12);

rad2deg=57.295779513082320876798154814105;
%-------偏航角解算为真航向---
psi_hmr=psi*rad2deg;
if (psi_hmr<   0.0) psi_hmr=psi_hmr+360.0; end
if (psi_hmr>=360.0) psi_hmr=psi_hmr-360.0; end
%-------输出量--------------
sys = [Vt;alpha*rad2deg;beta*rad2deg;phi*rad2deg;theta*rad2deg;psi*rad2deg;P*rad2deg;Q*rad2deg;R*rad2deg;PN;PE;H;(theta-alpha)*rad2deg];
% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

% 求导
% ---- 计算推力----
Pow =eng/100*(mass*g/3.0);
% ---- 计算空气密度和动压----
[ru,mach] = UAV_density(H,Vt);   % [air density] [mach number]
qs = SA*(ru*Vt*Vt/2);         % [Dynamic pressure](kg/m^2)
% ---- 计算气动力----
[CL0,CL_ele]=UAV_CL(alpha_deg);
L = qs*(CL0 + CL_ele*ele);
[CD0]=UAV_CD(alpha_deg);
D = qs*(CD0);

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate

%=============================================================================
%UAV_density密度和马赫数计算
%air density and mach number respect to the altitude and airspeed
%=============================================================================
function [ru,mach]=UAV_density(H,Vt)

if H<11000
   temp=1-0.0225569*H/1000;
   ru  =0.12492*9.8*exp(4.255277*log(temp));
   mach=Vt/340.375/sqrt(temp);
else
   temp=11-H/1000;
   ru  =0.03718*9.8*exp(temp/6.318);
   mach=Vt/295.188;
end
% end UAV_density

%=============================================================================
%UAV_CL  气动升力系数/导数计算
%=============================================================================
function [CL0,CL_ele]=UAV_CL(alpha_deg)
    IDX_alpha =[-4;-2;0;  2;4;8;12;16;20];
    TBL_CL0 =[-0.219;-0.04;0.139;  0.299;0.455;0.766;1.083;1.409;1.743];
    CL0 =interp1d(TBL_CL0,IDX_alpha,alpha_deg);
    CL_ele =0.00636;
%end UAV_CL
%=============================================================================
%UAV_CD  气动阻力系数/导数计算
%=============================================================================
function [CD]=UAV_CD(alpha_deg)
    IDX_alpha =[-4;-2;0;  2;4;8;12;16;20];
    TBL_CD0 =[ 0.026;0.024;0.024;0.028;0.036;0.061;0.102;0.141;0.173];
    CD =interp1d(TBL_CD0,IDX_alpha,alpha_deg);
%end UAV_CD
%=============================================================================
%interp1d  一维插值
%=============================================================================
function  y = interp1d(A,idx,xi)

if xi<idx(1)
   r=1;
elseif xi<idx(end)
   r = max(find(idx <= xi));      
else
   r = length(idx)-1;
end

DA = (xi-idx(r))/(idx(r+1)-idx(r));
y = A(r) + (A(r+1)-A(r))*DA;
% END interp1d

%
function [CM0, CM_alpha] = UAV_CM(alpha_deg)
    IDX_alpha =[-4; -2; 0; 2; 4; 8; 12; 16; 20];
    TBL_CM0 =[0.1161; 0.0777; 0.0393; 0.0009; -0.0375; -0.0759; -0.1527; -0.2295; -0.3063];
    CM0 =interp1d(TBL_CM0,IDX_alpha,alpha_deg);
    CM_alpha = -0.0192;
%

