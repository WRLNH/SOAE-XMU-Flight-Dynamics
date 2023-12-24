function [CL0, CL_ele] = UAV_CL(alpha_deg)
    % 此处传入的角度单位为角度制

    IDX_alpha = [-4; -2; 0; 2; 4; 8; 12; 16; 20];
    TBL_CL0 = [-0.219; -0.04; 0.139; 0.299; 0.455; 0.766; 1.083; 1.409; 1.743];
    CL0 = interp1d(TBL_CL0, IDX_alpha, alpha_deg);
    CL_ele = 0.00636;
end
