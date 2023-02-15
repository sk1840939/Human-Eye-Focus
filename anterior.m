function [results] = anterior(file, degree)
    %% Read & Process Data
    % Read data
    data_path = strcat("data/", file);

    % Reads all X/Y data from correct sheet and stores in matrix M. 
    % First column = x, Second column = y
    M = readmatrix(data_path, 'Sheet', 'Centered and Aligned', 'Range', 'A:B');
    figure; scatter(M(:,2), M(:,1)); title("Original plot of data"); xlim([-5,5]); ylim([-5,5]); % Plot to confirm

    % Bottom in original data, is indented from suture, so I
    % replicate from top side to smoothen out data and make consistent
    ant = M(M(:, 1) < 0, :); % Filter out anterior
    figure; scatter(ant(:,1), ant(:,2)); title("Anterior Raw"); % Plot to check

    % Process data for algorithm (places anterior on top, optic axis is x-axis)
    ant_new = sortrows(ant, 2);
    X_data = ant_new(:,2); Y_data = -1*ant_new(:,1);

    figure; plot(X_data, Y_data, 'LineWidth', 2); hold on;

    %% Fourier 
    fourier_bounds = [pi/2, 3*pi/2]; % Polar coordinates

    % Find fit of fourier
    [x_fourier, y_fourier] = fourier_fit(M(:,1), M(:,2));
    temp = x_fourier; x_fourier = y_fourier; y_fourier = -temp; % switching x & y to properly orient anterior surface
    fp = fplot(x_fourier, y_fourier, fourier_bounds, 'LineWidth', 2); X_fourierAnt = fp.XData; Y_fourierAnt = fp.YData;

    a_ant_fourier = max(X_fourierAnt);

    %% Chien - fit to raw
    syms t
    
    % Filter data
    chien_data = [X_data, Y_data];
%     bound = max(X_data); 
    
%     chien_data = chien_data(-bound < chien_data(:,1), :);
%     chien_data = chien_data(chien_data(:,1) < bound, :);
    offset = min(chien_data(:,2));
    X_chien = chien_data(:, 1); Y_chien = chien_data(:, 2) - offset; % Subtract offset for curvefitting
    
    % Used for ellipse fitting
    a_chien = max(X_chien) + 0.0001; % add epsilon for numerical stability
    a_ant = max(X_data)+ 0.0001;
    b0_ant = max(Y_data);

    % Get chien equation
    [x_chienAnt, y_chienAnt] = ChienAnt_CurveFive(X_chien, Y_chien);

    y_chienAnt = y_chienAnt + offset; % Re-add offset

    chien_bounds = [-pi/2, pi/2];
    fp = fplot(x_chienAnt, y_chienAnt, chien_bounds, 'LineWidth', 2); X_chienAnt = fp.XData; Y_chienAnt = fp.YData;

    %% Forbes
    % format data to forbes specs
    Y_forbes = -1*Y_data + max(Y_data); %figure; scatter(X_data, Y_forbes);

    % fit forbes to data
    syms rho;
    [forbes_eq, Y_forbes_raw, A, vertex_curvature] = forbes_severe(X_data', Y_forbes', degree);
    vroc = 1/vertex_curvature;
    
    forbes_reformat = -1*forbes_eq + double(vpa(subs(forbes_eq, rho, a_ant))); % Invert and offset
    forbes_eq = forbes_reformat + abs(max(Y_data) - max(Y_forbes_raw)); % Offset

    fp = fplot(rho, forbes_eq, [min(X_data), max(X_data)], 'LineWidth', 2); X_forbes = fp.XData; Y_forbes = fp.YData;

    forbes_eq = subs(forbes_eq, rho, t); % Change of variables for uniformity
    
    % Note -- the t in the forbes equation is cartesian! Stands for x (not
    % theta)
    
    %% Ellipse
    x_elipAnt = a_ant*cos(t); % in mm
    y_elipAnt = b0_ant*sin(t); % in mm

    elip_bounds = [0 pi];
    fp = fplot(x_elipAnt, y_elipAnt, elip_bounds, 'LineWidth', 2); X_elipAnt = fp.XData; Y_elipAnt = fp.YData;

    set(gca,'TickDir','out'); ax=gca; ax.FontSize=16;
    legend("Raw", "Fourier", "Chien", "Forbes", "Ellipse", 'Location', 'southwest', 'FontSize', 20); title("Raw Data and the 4 Analytical Models", 'FontSize', 24); xlabel("(mm)", 'FontSize', 20); ylabel("(mm)", 'FontSize', 20);

    %% Metrics
    zone = 3; % optical zone [-zone, +zone]
    offset_chien = abs(chien_bounds(1)) - asin(zone/a_ant); % polar - difference to come in from edges
    offset_elip = -1* (abs(elip_bounds(1)) - acos(zone/a_ant));
    offset_fourier = abs(fourier_bounds(1)) - asin(zone/a_ant_fourier); 

    % Find curvature
    k_chien = findCurvature(x_chienAnt, y_chienAnt, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
    k_elip = findCurvature(x_elipAnt, y_elipAnt, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
    k_fourier = findCurvature(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
    k_forbes = findCurvature(t, forbes_eq, -zone, zone)

    % Plot curvature - convert x-axis from radians to cartesian
    figure; hold on;
    fp = fplot(a_ant*sin(t), abs(k_chien), chien_bounds, 'LineWidth', 1); K_chien_x = fp.XData; K_chien_y = fp.YData;
    fp = fplot(a_ant*sin(t-pi/2), abs(k_elip), elip_bounds, 'LineWidth', 1); K_elip_x = fp.XData; K_elip_y = fp.YData;
    fp = fplot(a_ant_fourier*sin(t-pi), abs(k_fourier), fourier_bounds, 'LineWidth', 1); K_fourier_x = fp.XData; K_fourier_y = fp.YData;
    fp = fplot(t, abs(k_forbes), 'LineWidth', 1); K_forbes_x = fp.XData; K_forbes_y = fp.YData;
    set(gca,'TickDir','out'); ax=gca; ax.FontSize=16;
    legend("Chien", "Ellipse", "Fourier", "Forbes", 'FontSize', 20); title("Anterior Curvature for the 20-Year Old Lens", 'FontSize', 24); xlabel("Optical Zone Radius (mm)", 'FontSize', 20); ylabel("Magnitude of Curvature (1/mm)", 'FontSize', 20); ylim([0, 0.8]); xlim([-3,3]);

    % Find smoothing energy (integral of derivative of curvature squared)
    smth_chienAnt = double(vpa(vpaintegral(diff(k_chien, t, 1) ^ 2, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien)));
    smth_elipAnt = double(vpa(vpaintegral(diff(k_elip, t, 1) ^ 2, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip)));
    smth_fourierAnt = double(vpa(vpaintegral(diff(k_fourier, t, 1) ^ 2, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier)));
    smth_forbes = double(vpa(vpaintegral(diff(k_forbes, t, 1) ^ 2, -zone, zone)))

    % Find bending energy
    [bendE_chienAnt, firstD_chienAnt, expr_chienAnt] = findBendingEnergy(x_chienAnt, y_chienAnt, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
    [bendE_elipAnt, firstD_elipAnt, expr_elipAnt] = findBendingEnergy(x_elipAnt, y_elipAnt, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
    [bendE_fourierAnt, firstD_fourierAnt, expr_fourierAnt] = findBendingEnergy(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
    [bendE_forbes, firstD_forbes, expr_forbes] = findBendingEnergy(t, forbes_eq, -zone, zone)

    % Arc Length (over entire surface)
    length_chienAnt = double(vpaintegral(sqrt(diff(x_chienAnt, t, 1)^2 + diff(y_chienAnt, t, 1)^2), chien_bounds(1), chien_bounds(2)));
    length_elipAnt = double(vpaintegral(sqrt(diff(x_elipAnt, t, 1)^2 + diff(y_elipAnt, t, 1)^2), elip_bounds(1), elip_bounds(2)));
    length_fourierAnt = double(vpaintegral(sqrt(diff(x_fourier, t, 1)^2 + diff(y_fourier, t, 1)^2), fourier_bounds(1), fourier_bounds(2)));
    length_forbes = double(vpaintegral(sqrt(diff(t, t, 1)^2 + diff(forbes_eq, t, 1)^2), min(X_data), max(X_data)))
    
    % Calculate next set of metrics for the following optical zones 
    zone_vec = [0.5, 1, 1.5, 2, 2.5, 3]
    valROC_results = zeros(6, 4);
    varROC_results = zeros(1,4);
    for i = [1:6]
        i_zone = zone_vec(i);
        
        offset_chien = abs(chien_bounds(1)) - asin(i_zone/a_chien); % polar - difference to come in from edges
        offset_elip = -1* (abs(elip_bounds(1)) - acos(i_zone/a_ant));
        offset_fourier = abs(fourier_bounds(1)) - asin(i_zone/a_ant_fourier); 
        
        % Variance & Values of RoC
        n = linspace(0, chien_bounds(2) - offset_chien);
        yROC_chien = abs(double(vpa(subs(k_chien, n))));
        varROC_chien = var(yROC_chien);
        valROC_chien = abs(double(vpa(subs(1/k_chien, t, chien_bounds(1)+offset_chien))))

        n = linspace(0, elip_bounds(2) - offset_elip);
        yROC_elip = abs(double(vpa(subs(k_elip, n))));
        varROC_elip = var(yROC_elip);
        valROC_elip = abs(double(vpa(subs(1/k_elip, t, elip_bounds(1)+offset_elip))))

        n = linspace(0, fourier_bounds(2)-offset_fourier);
        yROC_fourier = abs(double(vpa(subs(k_fourier, n))));
        varROC_fourier = var(yROC_fourier);
        valROC_fourier = abs(double(vpa(subs(1/k_fourier, t, fourier_bounds(1)+offset_fourier))))
   
        n = linspace(0, i_zone);
        yROC_forbes = abs(double(vpa(subs(k_forbes, n))));
        varROC_forbes = var(yROC_forbes);
        valROC_forbes = abs(double(vpa(subs(1/k_forbes, t, -i_zone))))
        
        % Save value of ROC for each zone
        valROC_results(i, :) = [valROC_chien, valROC_forbes, valROC_fourier, valROC_elip];
        
        % Only save variance of ROC for zone is 3
        if i_zone == 3
            varROC_results = [varROC_chien, varROC_forbes, varROC_fourier, varROC_elip];
        end
    end

    % Fit (in microns)
    data = [X_data, Y_data];
    data_fit = data(-zone < data(:,1), :);
    data_fit = data_fit(data_fit(:,1) < zone, :);
    X_fit = data_fit(:, 1); Y_fit = data_fit(:, 2);

    % Converts from cartesian to radians 
    fit_forbes = getFit(X_fit, Y_fit, forbes_eq);
    fit_elip = getFit(asin(X_fit ./ a_ant)+pi/2, Y_fit, y_elipAnt);
    fit_chien = getFit(asin(X_fit ./ a_chien), Y_fit, y_chienAnt);
    fit_fourier = getFit(atan2(Y_fit, X_fit)+pi/2, Y_fit, y_fourier);
    
    metric_results = [fit_chien, fit_forbes, fit_fourier, fit_elip; 
                      length_chienAnt, length_forbes, length_fourierAnt, length_elipAnt;
                      vertex_curvature, vertex_curvature, vertex_curvature, vertex_curvature;
                      bendE_chienAnt, bendE_forbes, bendE_fourierAnt, bendE_elipAnt;
                      smth_chienAnt, smth_forbes, smth_fourierAnt, smth_elipAnt;
                      varROC_results];
    
    results = [metric_results; valROC_results];  
    
    % Plot first derivative of curvature
%     figure; hold on;
%     fplot(a_ant*sin(t), diff(k_chien)^2, chien_bounds, 'LineWidth', 1);
%     fplot(a_ant*sin(t-pi/2), diff(k_elip)^2, elip_bounds, 'LineWidth', 1);
%     fplot(a_ant_fourier*sin(t-pi), diff(k_fourier)^2, fourier_bounds, 'LineWidth', 1);
%     fplot(t, diff(k_forbes)^2, 'LineWidth', 1);
%     set(gca,'TickDir','out'); ax=gca; ax.FontSize=16;
%     legend("Chien", "Ellipse", "Fourier", "Forbes", 'FontSize', 20); title("Square of the First Derivative of the Anterior Curvature for the 20-Year Old Lens", 'FontSize', 24); xlabel("Optic Zone Radius (mm)", 'FontSize', 20); ylabel("First Derivative of Curvature Squared (1/mm^4)", 'FontSize', 20); ylim([0, 0.8]); xlim([-3,3]);
%     figure
    
%     curve_title = ["Raw X", "Raw Y", "Chien X", "Chien Y", "Ellipse X", "Ellipse Y", "Fourier X", "Fourier Y", "Forbes X", "Forbes Y"];
%     curve_data = [X_data'; Y_data'; X_chienAnt,NaN(1,1121); Y_chienAnt,NaN(1,1121); X_elipAnt,NaN(1,1114); Y_elipAnt,NaN(1,1114); X_fourierAnt,NaN(1,1085); Y_fourierAnt,NaN(1,1085); X_forbes,NaN(1,1200); Y_forbes,NaN(1,1200)];
%     k_title = ["Curvature Chien X", "Curvature Chien Y", "Curvature Ellipse X", "Curvature Ellipse Y", "Curvature Fourier X", "Curvature Fourier Y", "Curvature Forbes X", "Curvature Forbes Y"];
%     k_data = [K_chien_x,NaN(1,355); K_chien_y,NaN(1,355); K_elip_x,NaN(1,402); K_elip_y,NaN(1,402); K_fourier_x; K_fourier_y; K_forbes_x,NaN(1,520); K_forbes_y,NaN(1,520)];
%         
%     writematrix(curve_title, "ant_graph_data.xls", "Sheet", "Anterior Raw Data", "Range", 'A1')
%     writematrix(curve_data', "ant_graph_data.xls", "Sheet", "Anterior Raw Data", "Range", 'A2')
%     
%     writematrix(k_title, "ant_graph_data.xls", "Sheet", "Anterior Curvature Data", "Range", 'A1')
%     writematrix(k_data', "ant_graph_data.xls", "Sheet", "Anterior Curvature Data", "Range", 'A2')

    
end