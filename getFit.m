function [fit] = getFit(raw_x, raw_y, eq)
% Finds fit of equation to raw data
syms t
syms y(t)

y(t) = eq;

eq_y = double(y(raw_x)); % Value of equation evaluated at all points

% figure; hold on;
% scatter(raw_x, eq_y)
% scatter(raw_x, raw_y)
% legend('curve', 'raw')

% Calculate root mean-sqaured fit (and convert from mm to microns)
fit = sqrt(mean((eq_y - raw_y) .^ 2)) * 10^3;
end

