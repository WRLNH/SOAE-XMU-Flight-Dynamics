function [CY] = UAV_CY(beta_deg)
    CY_beta = -0.00909;
    % rad2deg = 57.295779513082320876798154814105;
    % beta_deg = beta * rad2deg;
    CY = beta_deg * CY_beta;
end
