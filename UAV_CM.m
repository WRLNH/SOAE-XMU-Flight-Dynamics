function [CM_alpha, CM_ele, CM_Q, CM_dalpha] = UAV_CM(alpha_deg)
    % 此处传入的角度单位为角度制

    IDX_alpha = [-4; -2; 0; 2; 4; 8; 12; 16; 20];
    TBL_CM0 = [0.1161; 0.0777; 0.0393; 0.0009; -0.0375; -0.0759; -0.1527; -0.2295; -0.3063];
    CM_alpha = interp1d(TBL_CM0, IDX_alpha, alpha_deg);
    CM_ele = -0.02052;
    CM_Q = -9.3136;
    CM_dalpha = -0.0192;
end
