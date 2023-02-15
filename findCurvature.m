function [k] = findCurvature(x, y, a, b)

%% Finds curvature
syms t;
% pre-calculating first/second order derivatives of x(t), y(t)
x_p = diff(x, t, 1);
x_pp = diff(x, t, 2);
y_p = diff(y, t, 1);
y_pp = diff(y, t, 2);

% formula for a curvature of a curve parameterized by x(t), y(t) -
% can be found on wikipedia
k = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);

end