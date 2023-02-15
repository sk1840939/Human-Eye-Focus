function [x_chien, y_chien] = ChienAnt_CurveFive(X, Y)
% Finds coefficients of Curve 5 from Type 4 of chien using least-squares
    
    % Definitions
    a = max(X);
    b = max(Y);
    theta = asin(X./a);
    
    % Set up and solve least-squares
    A = [b.*ones(size(theta)).* cos(theta), theta.^2.* cos(theta), theta.^4.* cos(theta)];
    coef = ((A'*A)\A') * Y;
    
    b0 = coef(1); b1 = coef(2); b3 = coef(3);
    
    % Construct curve
    syms t
    y_chien = (b + b1*t^2 + b3*t^4) * cos(t);
    x_chien = a*sin(t);   
    
end