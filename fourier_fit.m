function [x_param, y_param] = fourier_fit(X_data, Y_data)
% Use least-sqaures to fit the fourier equation to data passed in    

    % Convert to polar
    theta = atan2(Y_data, X_data);
    rho = sqrt(X_data.^2 + Y_data.^2);
    
    % Define function
    fun = @(x, t)x(1)*cos(0) + x(2)*cos(t) + x(3)*cos(2*t) + x(4)*cos(3*t) + x(5)*cos(4*t) + x(6)*cos(5*t) + x(7)*cos(6*t) + x(8)*cos(7*t) + x(9)*cos(8*t) + x(10)*cos(9*t) + x(11)*cos(10*t);
    x0 = ones(1,11);
    
    % Determine Fourier equation coefficients
    coef = lsqcurvefit(fun, x0, theta, rho);
    
    syms t
    
    % Reconstruct fourier equation
    sym_eq = sym(fun(coef, t));
    
    % Convert back to cartesian
    x_param = sym_eq * cos(t);
    y_param = sym_eq * sin(t);
    
end