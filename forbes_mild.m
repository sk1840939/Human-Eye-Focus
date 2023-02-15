function [z, Y_forbes, A] = forbes_mild(X_data, Y_data, M)
    %% Define test data & input args
%     M = 32;
% 
%     syms u;
%     y = 1/40 * u^2;
%     fp = fplot(y, [-20, 20], 'MeshDensity', 100); X_data = fp.XData; Y_data = fp.YData;

    %% Construct p
    syms x
    % Initial conditions
    P = [2, 6-8*x];

    for i = 3:M+1
        P(i) = (2-4*x) * P(i-1) - P(i-2);
    end

    %% Determine g, h, and f
    % Initial conditions
    g = [-1/2];
    h = [];
    f = [2,sqrt(19)/2];

    % i used to index, m used for actual value
    for i = 3:M+3
        m = i-1;

        h(i-2) = (-m*(m-1)) / (2*f(i-2));

        if i <= M+2
            g(i-1) = -(1 + g(i-2)*h(i-2)) / f(i-1);
        end

        if i <= M+1
            f(i) = sqrt(m*(m+1) + 3 - g(i-1)^2 - h(i-2)^2);
        end
    end

    %% Determine Q
    % Initial conditions
    Q = [1, (13-16*x)/sqrt(19)];

    for i = 3:M+1
        Q(i) = (P(i) - g(i)*Q(i-1) - h(i-2)*Q(i-2)) / f(i);    
    end
    % It can be verified that these basis values match (14) from Axis.

    %% Solve for A 
    % Curvature of best-fit circle
    c = (2*max(Y_data) / (max(X_data)^2 + max(Y_data)^2)); 

    Y = Y_data'; % (dx1)

    % Normalize data
    rho = X_data;
    norm_data = rho ./ max(rho);

    % Generate Q with data
    Q_data = double(vpa(subs(Q, norm_data'.^2)))'; % (M x d) % M rows, d cols

    % Assemble feature matrix (S)
    lambda = (c * rho.^2) ./ (1+sqrt(1-c^2 .* rho.^2));
    beta = (norm_data.^2 .* (1 - norm_data.^2)) ./ (sqrt(1-c^2.*rho.^2));

    S = [lambda ; beta .* Q_data];

    A = ((S*S')\S) * Y;

    %% Assembling curve 

    % Points
    Y_forbes = A'*S;

    % Equation (rho is equivalent to x, z is equivalent to y)
    syms rho
    u = rho / max(X_data);

    z = (c*rho^2) / (1+sqrt(1-c^2*rho^2)) + (u^2 * (1-u^2)) / (sqrt(1-c^2*rho^2)) .* (A(2:end)' * subs(Q, u^2)');
end