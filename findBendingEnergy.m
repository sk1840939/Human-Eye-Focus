function [bend_E, first_d, expr] = findBendingEnergy(x_eq, y_eq, a, b)
%% Finds bending energy
syms t;

% First derivatives wrt t
y_p = diff(y_eq, t, 1);
x_p = diff(x_eq, t, 1);

% Numerator
num = diff(y_p/x_p, t, 1);

% Final Integral Expression
y_pp_x = num / x_p;
expr = y_pp_x ^ 2;

% Evaluate integral from bounds
bend_E = eval(vpaintegral(expr, a, b));

first_d = eval(vpaintegral((y_p/x_p) ^ 2, a, b));
end