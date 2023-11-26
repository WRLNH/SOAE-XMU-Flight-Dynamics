function [CD0] = UAV_CD(alpha_deg)
    % 此处传入的角度单位为角度制

    IDX_alpha = [-4; -2; 0; 2; 4; 8; 12; 16; 20];
    TBL_CD0 = [0.026; 0.024; 0.024; 0.028; 0.036; 0.061; 0.102; 0.141; 0.173];
    % rad2deg = 180 / pi;
    % alpha_deg = alpha * rad2deg;
    CD0 = interp1d(TBL_CD0, IDX_alpha, alpha_deg);

end
