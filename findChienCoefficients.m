function [b1,b3,b4,b5] = findChienCoefficients(X, Y, a, b0)
%FINDCHIENCOEFFICIENTS 
% Finds coefficients for chien using NNLS and the coordinates from the 
% Fourier curve

% Defining curve to fit
min(X)
max(X)
theta = strcat('asin(x/', num2str(a), ')');
cos_theta = strcat('cos(', theta, ')');

b0_term = strcat(num2str(b0),'*',cos_theta);
b1_term = strcat('b*', theta, '^2*',cos_theta);
b3_term = strcat('c*', theta, '^4*', cos_theta);
%b4_term = strcat('d*', theta, '^5*', cos_theta);
%b5_term = strcat('e*', theta, '^6*', cos_theta);

ft = fittype(strcat(b0_term, '+', b1_term, '+', b3_term), 'dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'b', 'c'})
options = fitoptions('Method', 'NonLinearLeastSquares', 'StartPoint', [-0.1, 0.1]);

% Fit using NNLS
myfit = fit(X', Y', ft, options)

% Extract coefficients
b1 = myfit.b;
b3 = myfit.c;
%b4 = myfit.d;
%b5 = myfit.e;
end

