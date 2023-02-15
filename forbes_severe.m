function [z, Y_forbes, A, c] = forbes_severe(X_data, Y_data, M)
    %% Fit ellipse
    syms t
    a = max(X_data);
    b = max(Y_data);
    
    x_elip = a*cos(t); % in mm
    y_elip = b*sin(t) + b;
    vertex = 3*pi/2;
    e = sqrt(1-(b^2/a^2)); % eccentricity
    
    %% Fit hyperbola -- unused
%     y_hyp = b*sec(t)-b;
%     x_hyp = a*tan(t);
%     vertex = 0;
%     e = sqrt(1+(b^2/a^2)); % eccentricity
    
%     figure; hold on;
%     scatter(X_data, Y_data);
%     fplot(x_elip, y_elip, [pi, 2*pi]);

    %% Find conic constant and radius of curvature
    k = -e^2; % conic constant
    epsilon = 1+k; % Forbes's conic parameter
    
    elip_k = findCurvature(x_elip, y_elip, 0, 2*pi);
    c = double(vpa(subs(elip_k, t, vertex))); % Curvature at vertex of curve
%     fplot(t, elip_k)
    
    %% Determine Q
    syms x
    Q = jacobiP([0:M], 0, 4, 2*x-1);
    
    %% Solve for A 
    Y = Y_data'; % (dx1)

    % Normalize data
    rho = X_data;
    norm_data = rho ./ max(rho);

    % Generate Q with data
    Q_data = double(vpa(subs(Q, norm_data'.^2)))'; % (M x d) % M rows, d cols

    % Assemble feature matrix (S)
    lambda = (c * rho.^2) ./ (1+sqrt(1 - (epsilon.*c^2.*rho.^2))); % spheric component
    beta = (norm_data.^2 .* (1 - norm_data.^2)) ./ (sqrt(1-c^2.*rho.^2.*norm_data.^2)); % aspheric weighting on Jacobi polynomial

    S = [lambda ; beta .* Q_data];
%     S = [lambda ; (norm_data.^4) .* Q_data];

    % Weighting described in Forbes Method
    weights = (norm_data.*(1-norm_data)).^(-1/2);
    weights(isinf(weights)) = max(weights(weights < Inf));
    W = diag(abs(weights));

    % Weighted Least Squares to find coefficients
    A = ((S*W*S')\S) *W* Y;
    
%     scatter(X_data, lambda) % visualize

    %% Assembling curve 
    % Invert WLS to solve for curve
    Y_forbes = A'*S;
%     scatter(X_data, Y_forbes); 

    % Assemble Forbes Equation (rho is equivalent to x, z is equivalent to y)
    syms rho
    u = rho / max(X_data);

    lambda = (c*rho^2) /  (1+sqrt(1 - (epsilon*c^2*rho^2)));
    beta = (u.^2 .* (1 - u.^2)) ./ (sqrt(1-c^2.*rho^2*u.^2));
    
    S = [lambda; (beta .* subs(Q, u^2))];
    offset = abs(max(Y_data) - max(Y_forbes));
    z = (A' * S);
    
%     fplot(rho, z);
%     legend("Raw", "Ellipse", "Forbes (Points)", "Forbes (Curve)")
%     title("Forbes Severe");
end