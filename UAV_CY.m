function [CY] = UAV_CY(beta_deg)
    % 此处传入的角度单位为角度制

    CY_beta = -0.00909;
    CY = beta_deg * CY_beta;
end
