function [x_param, y_param] = fourier_curve(age)
%FOURIER_CURVE 
% Given an age, returns the parameterized syms equation for 
% the curve (in cartesian coordinates)

    syms t;
    % Given coefficients (paper) % original is 2.6466 %
    A_n1 = [3 0.2246 -0.97938 0.010573 0.37993 -0.032321 -0.16846 0.027934 0.066522 -0.014232 -0.021375];
    A_n2 = [812.11 170.62 -297.37 -34.901 -26.276 1.6647 69.192 -9.5571 -42.251 1.7295 18.638] .* 1e-5;

    A = A_n1 + A_n2 .* age;
    
    sum(A)

    % Defining polar equation
    rho = t-t;

    for i=0:10
        rho = rho + A(i+1)*cos(i*t); % could convert cos to sin using cos(t) = sin(pi/2 - t)
    end

    % Parameterize the curve (from polar to cartesian)
    % THEY ARE SWAPPED TO SWITCH ANTERIOR/POSTERIOR INTO FUNCTION (Flips by 90 degrees)
    y_param = -rho * cos(t); % x = r*cos(theta)
    x_param = rho * sin(t); % y = r*sin(theta)
    
    
    ds = sqrt(rho^2 + diff(rho, t, 1)^2);
    
    rho
    pretty(rho)
    
          
    
end

