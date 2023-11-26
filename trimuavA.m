function [xtrim, ytrim, utrim, dxtrim] = trimuavA(V, H, path)
    % [xtrim,utrim]=cktrim(V,H,path)
    % Calculate trim date
    %         1      2          3     4      5   6   7   8   9   10   11   12
    % x = [Vt, alpha, beta, PN, PE, alt,  P,  Q, R, phi, theta, psi,]
    %       1    2     3     4      5     6    7   8  9   10   11   12
    % Y = [Vt, alpha, beta, PN, PE, alt,  P,  Q, R, phi, theta, psi,theta-alpha]
    %  u = [elevator, aileron, rudder, engine]

    x0 = [V; 3/57.3; 0; 0; 0; H; 0; 0; 0; 0; 3/57.3 + path / 57.3; 0];
    IX = [1; 7; 8; 9; 10; 12];

    u0 = [0; 0; 0; 0.3];
    IU = [2; 3];

    y0 = [V; 0; 0; 0; 0; H; 0; 0; 0; 0; 0; 0; path];
    IY = [1; 7; 8; 9; 10; 12; 13];

    dx0 = zeros(12, 1);
    IDX = [1; 2; 3; 7; 8; 9; 10; 11; 12];

    options = optimset('MaxFunEvals', 2e4, 'TolX', 1e-4, 'TolFun', 1e-4);

    [xtrim, utrim, ytrim, dxtrim] = trim('uavCtrl', x0, u0, y0, IX, IU, IY, dx0, IDX);

    [a, b, c, d] = linmod('uavCtrl', xtrim, utrim);
    % Vt=x(1);alpha=x(2);beta=x(3);    PN=x(4);   PE=x(5);    H=x(6);
    % P=x(7);      Q=x(8);   R=x(9);  phi=x(10);theta=x(11);psi=x(12); path=13
    [along, blong, clong, dlong] = ssselect(a, b, c, d, [1; 4], [1; 2; 8; 11; 6], [1; 2; 8; 11; 6]);
    [alate, blate, clate, dlate] = ssselect(a, b, c, d, [2; 3], [3; 7; 9; 10; 12; 5], [3; 7; 9; 10; 12; 5]);
    save('data', 'xtrim', 'utrim', 'ytrim', 'dxtrim', 'a', 'b', 'c', 'd', 'along', 'blong', 'clong', 'dlong', 'alate', 'blate', 'clate', 'dlate');
    % damp(along);
    % damp(alate);
    % [num,den]=ss2tf(along,blong,clong(4,:),dlong(4,:),1);
    % printsys(num,den,'s');
end