function [sys,x0,str,ts] = uavL(t,x,u,flag)
%---------the author is HLJ。
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
rad2deg=180/pi;

psi_hmr=0.0;
if (psi_hmr> 180)   
    psi_hmr = psi_hmr-360.0; end
if (psi_hmr<-180)   
    psi_hmr = psi_hmr+360.0; end

Vt = 25; alpha = 0; beta = 0 + 1; theta = 0; H = 1000;
x0  = [Vt; alpha/rad2deg; beta/rad2deg;0; 0;H; 0;0;0; 0; theta/rad2deg;psi_hmr/rad2deg];
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
Vt      = x(1);   PN = x(4);   P = x(7);   phi   = x(10);       
alpha   = x(2);   PE = x(5);   Q = x(8);   theta = x(11);        
beta    = x(3);    H = x(6);   R = x(9);   psi   = x(12);
rad2deg=180/pi;
alpha_deg = alpha * rad2deg; beta_deg = beta * rad2deg; theta_deg = theta * rad2deg;

% ---输入量-------
ele = u(1);      % elevator deflection angle (deg)
ail = u(2);      % aileron  deflection angle (deg)
rud = u(3);      % rudder   deflection angle (deg)
eng = u(4);   
% ---质量数据-------
mass=17; Ix=1.71; Iy=5.54; Iz=4.15; cbar=0.423; SA=1.3536;b=3.2;g=9.81;


% ---- 计算推力----
Pow=eng/100*(mass*g)/4;
% ---- 计算空气密度和动压----
[ru,mach] = UAV_density(H,Vt);   % [air density] [mach number]
qs = SA*(ru*Vt*Vt/2);         % [Dynamic pressure](kg/m^2)
% ---- 计算气动力----
[CL0,CL_ele]=UAV_CL(alpha_deg);
L = qs*(CL0 + CL_ele*ele);
[CD]=UAV_CD(alpha_deg);
D = qs*CD;
[CY_beta]=UAV_CY(beta_deg);
Y = qs * CY_beta * beta_deg;


% ---- 基于气动数据计算气动受力，并以十二个方程求解状态量的微分量----
%----机体坐标系受力
Fx=Pow-D*cos(alpha)*cos(beta)-Y*cos(alpha)*sin(beta)+L*sin(alpha);%------L前面是正是负？
Fy=-D*sin(beta)+Y*cos(beta);
Fz=-D*sin(alpha)*cos(beta)-Y*sin(alpha)*sin(beta)-L*cos(alpha);
%-----机体坐标系速度------
U = Vt*cos(alpha)* cos(beta);
V = Vt*sin(beta);
W = Vt*sin(alpha)* cos(beta);

%----加速度
dU=R*V-Q*W-g*sin(theta)+Fx/mass;
dV=-R*U+P*W+g*sin(phi)*cos(theta)+Fy/mass;
dW=Q*U-P*V+g*cos(phi)*cos(theta)+Fz/mass;

dVt=(U*dU+V*dV+W*dW)/Vt;
dbeta=(dV*Vt-V*dVt)/(Vt*Vt*cos(beta)); 
dalpha=(U*dW-W*dU)/(U*U+W*W);


% ---- 计算气动力矩----
[CR0,CR_R,CR_P,CR_rud,CR_ail]=UAV_CR(beta_deg);%-----滚转
CR=CR0+b/(2*Vt)*(CR_R*R+CR_P*P)+CR_ail*ail+CR_rud*rud;
Lbar=qs*b*CR;
[CM_alpha,CM_ele,CM_Q,CM_dalpha]=UAV_CM(alpha_deg);%-----俯仰
CM=CM_alpha+CM_ele*ele+cbar/(2*Vt)*(CM_Q*Q);
M=qs*cbar*CM;
[CN0,CN_R,CN_P,CN_rud,CN_ail]=UAV_CN(beta_deg);%-----偏航
CN =CN0+b/(2*Vt)*(CN_R*R+CN_P*P)+CN_rud*rud+CN_ail*ail;
N=qs*b*CN;

dphi=P+tan(theta)*(Q*sin(phi)+R*cos(phi));
dtheta=Q*cos(phi)-R*sin(phi);
dpsi=1/cos(theta)*(Q*sin(phi)+R*cos(phi));

dP=1/(Ix*Iz)*(Iz*Lbar+(Iy-Iz)*Iz*R*Q);
dQ=1/Iy*((Iz-Ix)*P*R+M);
dR=1/(Ix*Iz)*((Ix-Iy)*Ix*P*Q+N*Ix);

dPN=U*cos(psi)*cos(theta)+V*(-sin(psi)*cos(phi)+cos(psi)*sin(theta)*sin(phi))+W*(sin(psi)*sin(phi)+cos(psi)*sin(theta)*cos(phi));
dPE=U*cos(theta)*sin(psi)+V*(cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi))+W*(-sin(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi));%-----bug
dH=U*sin(theta)-V*sin(phi)*cos(theta)-W*cos(phi)*cos(theta);

% ---状态量的微分量-------
sys = [dVt;dalpha;dbeta;dPN;dPE;dH;dP;dQ;dR;dphi;dtheta;dpsi;];
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
sys = [];
%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
% ---状态量-------
Vt      = x(1);   PN = x(4);   P = x(7);   phi   = x(10);       
alpha   = x(2);   PE = x(5);   Q = x(8);   theta = x(11);        
beta    = x(3);    H = x(6);   R = x(9);   psi   = x(12);   

rad2deg=57.295779513082320876798154814105;
%-------偏航角解算为真航向---
psi_hmr=psi*rad2deg;
if (psi_hmr<0.0) psi_hmr=psi_hmr+360.0; end
if (psi_hmr>=360.0) psi_hmr=psi_hmr-360.0; end
%-------输出量--------------
sys = [Vt;alpha*rad2deg;beta*rad2deg;PN;PE;H;P*rad2deg;Q*rad2deg;R*rad2deg;phi*rad2deg;theta*rad2deg;psi*rad2deg;(theta-alpha)*rad2deg];
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
function [CL0,CL_ele]=UAV_CL(alpha)
    IDX_alpha =[-4;-2;0;2;4;8;12;16;20];
    TBL_CL0 =[-0.219;-0.04;0.139; 0.299;0.455;0.766;1.083;1.409;1.743];
    rad2deg=57.295779513082320876798154814105;
    alpha_deg=alpha*rad2deg;
    CL0 =interp1d(TBL_CL0,IDX_alpha,alpha_deg);
    CL_ele =0.00636;
%end UAV_CL
%=============================================================================
%UAV_CD  气动阻力系数/导数计算
%=============================================================================
function [CD0]=UAV_CD(alpha)
    IDX_alpha =[-4;-2;0;  2;4;8;12;16;20];
    TBL_CD0 =[ 0.026;0.024;0.024;0.028;0.036;0.061;0.102;0.141;0.173];
    rad2deg=57.295779513082320876798154814105;
    alpha_deg=alpha*rad2deg;
    CD0 =interp1d(TBL_CD0,IDX_alpha,alpha_deg);
%end UAV_CD
%=============================================================================
%UAV_CY  气动侧力系数/导数计算
%=============================================================================
function [CY]=UAV_CY(beta)
   CY_beta=-0.00909;
   rad2deg=57.295779513082320876798154814105;
   beta_deg=beta*rad2deg;
   CY =beta_deg*CY_beta;
%end UAV_CY
%=============================================================================
%UAV_CR  滚转系数/导数计算
%=============================================================================
function [CR0, CR_R, CR_P,CR_rud,CR_ail]=UAV_CR(beta)
    CR_beta=-0.00600;
    CR_R=0.01832;
    CR_P=-0.52568;
    rad2deg=180/pi;
    beta_deg=beta*rad2deg;
    CR0 =CR_beta*beta_deg;
    CR_rud =-0.000144;
    CR_ail=-0.003618;
%end UAV_CR
%=============================================================================
%UAV_CM  俯仰系数/导数计算
%=============================================================================
function [CM_alpha,CM_ele,CM_Q,CM_dalpha]=UAV_CM(alpha)
    IDX_alpha =[-4;-2;0;2;4;8;12;16;20];
    TBL_CM0 =[0.1161;0.0777;0.0393; 0.0009;-0.0375;-0.0759;-0.1527;-0.2295;-0.3063];
    rad2deg=57.295779513082320876798154814105;
    alpha_deg=alpha*rad2deg;
    CM_alpha =interp1d(TBL_CM0,IDX_alpha,alpha_deg);
    CM_ele=-0.02052;
    CM_Q=-9.3136;
    CM_dalpha=-1.0354;
%end UAV_CM
%===================
% %UAV_CN  偏航系数/导数计算
%=============================================================================
function [CN0,CN_R,CN_P,CN_rud,CN_ail]=UAV_CN(beta)
    CN_beta=0.00235;
    CN_R=-0.15844;
    CN_P=-0.01832;
    rad2deg=180/pi;
    beta_deg=beta*rad2deg;
    CN0 =CN_beta*beta_deg;
    CN_rud =-0.00111;
    CN_ail =0.000132;
%end UAV_CN=
% =========================================================
%interp1d  一维插值
%=============================================================================
function  y = interp1d(A,idx,xi)

if xi<idx(1)
   r=1;
elseif xi<idx(end)
   r = max(find(idx <= xi));      
else
   r = length(idx)-1;%从零开始
end

DA = (xi-idx(r))/(idx(r+1)-idx(r));
y = A(r) + (A(r+1)-A(r))*DA;%一阶线性处理
% END interp1d



