
function  [results] = posterior(file, degree)
    %% Read & Process Data
    % Read data
    data_path = strcat("data/", file);

    % Reads all X/Y data from correct sheet and stores in matrix M. 
    % First column = x, Second column = y

    M = readmatrix(data_path, 'Sheet', 'Centered and Aligned', 'Range', 'A:B');
    figure; scatter(M(:,1), M(:,2)); title("Original plot of data"); xlim([-5,5]); ylim([-5,5]); % Plot to confirm

    % Bottom in original data, is indented from suture, so I
    % replicate from top side to smoothen out data and make consistent
    post = M(M(:, 1) > 0, :); % Filter out posterior
    figure; scatter(post(:,1), post(:,2)); title("Posterior Raw"); % Plot to check

%     post_top = post(post(:,2) > 0, :); % Filter out top of posterior
%     post_bot = [post_top(:,1), -1*post_top(:,2)]; % Flip across x-axis
% 
%     post_new = cat(1, post_top, post_bot); % Concat to form new posterior
%     figure; scatter(post_new(:,1), post_new(:,2)); title("Posterior w/ Fixed Suture") % Plot to check
%     
%     post_new = sortrows(post_new, 2);
%     X = post_new(:,1); Y = post_new(:,2); % Draw out X, Y from posterior data

    % Process data for algorithm (places posterior on top, optic axis is x-axis)
    post_new = sortrows(post, 2);
    X_data = post_new(:,2); Y_data = post_new(:,1);
    figure; plot(X_data, Y_data, 'LineWidth', 2); hold on;

    %% Fourier 
    % Replicate top of lens to remove suture so we can fit the fourier data
%     data_top = M(M(:,2) > 0, :);
%     data_bottom = [data_top(:,1), -1*data_top(:,2)];
% 
%     fourier_data = cat(1, data_top, data_bottom);
    
    % Find fit of fourier
    [x_fourier, y_fourier] = fourier_fit(M(:,1), M(:,2));
    temp = x_fourier; x_fourier = y_fourier; y_fourier = temp;

    % Obtain x,y coords for posterior
    fourier_bounds = [-pi/2, pi/2];
    fp = fplot(x_fourier, y_fourier, fourier_bounds); X_fourier = fp.XData; Y_fourier = fp.YData;

    a_post_fourier = max(X_fourier);

    %% Chien - fit to raw
    syms t
    
    chien_data = [X_data, Y_data];
    bound = max(X_data);
%     bound = 3.25; %max(X_data);
    chien_data = chien_data(-bound < chien_data(:,1), :);
    chien_data = chien_data(chien_data(:,1) < bound, :);
    offset = min(chien_data(:,2));
    X_chien = chien_data(:, 1); Y_chien = chien_data(:, 2) - offset;
    
    b0_chien = max(Y_chien);
    a_chien = max(X_chien) + 0.0001; % add epsilon for numerical stability
    a_post = max(X_data)+ 0.0001;
    b0_post = max(Y_data);

    [x_chien, y_chien] = ChienAnt_CurveFive(X_chien, Y_chien);

    y_chien = y_chien + offset;

    chien_bounds = [-pi/2, pi/2];
    fp = fplot(x_chien, y_chien, chien_bounds, 'LineWidth', 2); X_chien = fp.XData; Y_chien = fp.YData;

    %% Forbes
    % format data to forbes specs
    Y_forbes = -1*Y_data + max(Y_data); %figure; scatter(X_data, Y_forbes);

    % fit forbes to data
    syms rho;
    [forbes_eq, Y_forbes_raw, A, vertex_curvature] = forbes_severe(X_data', Y_forbes', degree);
    vroc = 1/vertex_curvature;    
    
    forbes_reformat = -1*forbes_eq + max(Y_data);
    forbes_eq = forbes_reformat

    fp = fplot(rho, forbes_eq, [min(X_data), max(X_data)]); X_forbes = fp.XData; Y_forbes = fp.YData;

    forbes_eq = subs(forbes_eq, rho, t);
    % Note -- the t in the forbes equation is cartesian! Stands for x (not
    % theta)
    %% Ellipse
    x_elip = a_post*cos(t); % in mm
    y_elip = b0_post*sin(t); % in mm

    elip_bounds = [0 pi];
    fp = fplot(x_elip, y_elip, elip_bounds); X_elip = fp.XData; Y_elip = fp.YData;
    %plot(X_elip, Y_elip, 'LineWidth', 2);

    set(gca,'TickDir','out'); ax=gca; ax.FontSize=16;
    legend("Raw", "Fourier", "Chien", "Forbes", "Ellipse", 'Location', 'southwest', 'FontSize', 20); title("Raw Data and the 4 Analytical Models", 'FontSize', 24); xlabel("(mm)", 'FontSize', 20); ylabel("(mm)", 'FontSize', 20);

    %% Metrics
    zone = 3; % optical zone [-zone, +zone]
    offset_chien = abs(chien_bounds(1)) - asin(zone/a_post); % polar - difference to come in from edges
    offset_elip = -1* (abs(elip_bounds(1)) - acos(zone/a_post));
    offset_fourier = abs(fourier_bounds(1)) - asin(zone/a_post_fourier); 

    % Find curvature
    k_chien = findCurvature(x_chien, y_chien, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
    k_elip = findCurvature(x_elip, y_elip, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
    k_fourier = findCurvature(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
    k_forbes = findCurvature(t, forbes_eq, -zone, zone)

    % Plot curvature - convert x-axis from radians to cartesian
    figure; hold on;
    fp = fplot(a_post*sin(t), abs(k_chien), chien_bounds); K_chien_x = fp.XData; K_chien_y = fp.YData;
    fp = fplot(a_post*sin(t-pi/2), abs(k_elip), elip_bounds); K_elip_x = fp.XData; K_elip_y = fp.YData;
    fp = fplot(a_post_fourier*sin(t-pi), abs(k_fourier), fourier_bounds); K_fourier_x = fp.XData; K_fourier_y = fp.YData;
    fp = fplot(t, abs(k_forbes)); K_forbes_x = fp.XData; K_forbes_y = fp.YData;
    set(gca,'TickDir','out'); ax=gca; ax.FontSize=16;
    legend("Chien", "Ellipse", "Fourier", "Forbes", 'FontSize', 20); title("Posterior Curvature for the 20-Year Old Lens", 'FontSize', 24); xlabel("Optical Zone Radius (mm)", 'FontSize', 20); ylabel("Magnitude of Curvature (1/mm)", 'FontSize', 20); ylim([0, 0.8]); xlim([-3,3]);
    figure

    % Find smoothing energy (integral of derivative of curvature squared)
    smth_chien = double(vpa(vpaintegral(diff(k_chien, t, 1) ^ 2, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien)));
    smth_elip = double(vpa(vpaintegral(diff(k_elip, t, 1) ^ 2, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip)));
    smth_fourier = double(vpa(vpaintegral(diff(k_fourier, t, 1) ^ 2, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier)));
    smth_forbes = double(vpa(vpaintegral(diff(k_forbes, t, 1) ^ 2, -zone, zone)))

    % Find bending energy
    [bendE_chien, firstD_chienAnt, expr_chienAnt] = findBendingEnergy(x_chien, y_chien, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien);
    [bendE_elip, firstD_elipAnt, expr_elipAnt] = findBendingEnergy(x_elip, y_elip, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip);
    [bendE_fourier, firstD_fourierAnt, expr_fourierAnt] = findBendingEnergy(x_fourier, y_fourier, fourier_bounds(1)+offset_fourier, fourier_bounds(2)-offset_fourier);
    [bendE_forbes, firstD_forbes, expr_forbes] = findBendingEnergy(t, forbes_eq, -zone, zone)

    % Arc Length (over entire surface)
    length_chien = double(vpaintegral(sqrt(diff(x_chien, t, 1)^2 + diff(y_chien, t, 1)^2), chien_bounds(1), chien_bounds(2)));
    length_elip = double(vpaintegral(sqrt(diff(x_elip, t, 1)^2 + diff(y_elip, t, 1)^2), elip_bounds(1), elip_bounds(2)));
    length_fourier = double(vpaintegral(sqrt(diff(x_fourier, t, 1)^2 + diff(y_fourier, t, 1)^2), fourier_bounds(1), fourier_bounds(2)));
    length_forbes = double(vpaintegral(sqrt(diff(t, t, 1)^2 + diff(forbes_eq, t, 1)^2), min(X_data), max(X_data)))
    
    zone_vec = [0.5, 1, 1.5, 2, 2.5, 3]
    valROC_results = zeros(6, 4);
    varROC_results = zeros(1,4);
    for i = [1:6]
        i_zone = zone_vec(i);
        
        offset_chien = abs(chien_bounds(1)) - asin(i_zone/a_chien); % polar - difference to come in from edges
        offset_elip = -1* (abs(elip_bounds(1)) - acos(i_zone/a_post));
        offset_fourier = abs(fourier_bounds(1)) - asin(i_zone/a_post_fourier); 
        
        % Mean/Variance of RoC - variance found numerically
%         meanROC_chien = abs(1/((chien_bounds(2)-offset_chien) - (chien_bounds(1)+offset_chien)) * vpa(vpaintegral(k_chien, chien_bounds(1)+offset_chien, chien_bounds(2)-offset_chien)));
        n = linspace(0, chien_bounds(2) - offset_chien);
        yROC_chien = abs(double(vpa(subs(k_chien, n))));
        varROC_chien = var(yROC_chien);
        valROC_chien = abs(double(vpa(subs(1/k_chien, t, chien_bounds(1)+offset_chien))))

%         meanROC_elip = abs(1/((elip_bounds(2)-offset_elip) - (elip_bounds(1)+offset_elip)) * vpa(vpaintegral(k_elip, elip_bounds(1)+offset_elip, elip_bounds(2)-offset_elip)));
        n = linspace(0, elip_bounds(2) - offset_elip);
        yROC_elip = abs(double(vpa(subs(k_elip, n))));
        varROC_elip = var(yROC_elip);
        valROC_elip = abs(double(vpa(subs(1/k_elip, t, elip_bounds(1)+offset_elip))))

%         meanROC_fourier = abs(1/((fourier_bounds(2)-offset_fourier) - (fourier_bounds(1)+offset_fourier)) * int_k);
        n = linspace(0, fourier_bounds(2)-offset_fourier);
        yROC_fourier = abs(double(vpa(subs(k_fourier, n))));
        varROC_fourier = var(yROC_fourier);
        valROC_fourier = abs(double(vpa(subs(1/k_fourier, t, fourier_bounds(1)+offset_fourier))))
   
        %meanROC_forbesAnt = abs(1/(2*zone) * vpa(vpaintegral(k_forbes, -zone, zone)));
        n = linspace(0, i_zone);
        yROC_forbes = abs(double(vpa(subs(k_forbes, n))));
        varROC_forbes = var(yROC_forbes);
        valROC_forbes = abs(double(vpa(subs(1/k_forbes, t, -i_zone))))
        
        valROC_results(i, :) = [valROC_chien, valROC_forbes, valROC_fourier, valROC_elip];
        if i_zone == 3
            varROC_results = [varROC_chien, varROC_forbes, varROC_fourier, varROC_elip];
        end
    end

    % Fit (in microns)
    data = [X_data, Y_data];
    data_fit = data(-zone < data(:,1), :);
    data_fit = data_fit(data_fit(:,1) < zone, :);
    X_fit = data_fit(:, 1); Y_fit = data_fit(:, 2);
    figure; hold on; scatter(X_fit, Y_fit); title("Data to fit");

    % Converts from cartesian to radians 
    fit_forbes = getFit(X_fit, Y_fit, forbes_eq);
    fit_elip = getFit(asin(X_fit ./ a_post)+pi/2, Y_fit, y_elip);
    fit_chien = getFit(asin(X_fit ./ a_chien), Y_fit, y_chien);
    fit_fourier = getFit(atan2(Y_fit, X_fit)+3*pi/2, Y_fit, y_fourier);
    
    metric_results = [fit_chien, fit_forbes, fit_fourier, fit_elip; 
                      length_chien, length_forbes, length_fourier, length_elip;
                      vertex_curvature, vertex_curvature, vertex_curvature, vertex_curvature;
                      bendE_chien, bendE_forbes, bendE_fourier, bendE_elip;
                      smth_chien, smth_forbes, smth_fourier, smth_elip;
                      varROC_results];
    
    results = [metric_results; valROC_results]; 
    
    % Plot first derivative of curvature
    figure; hold on;
    fplot(a_post*sin(t), diff(k_chien)^2, chien_bounds, 'LineWidth', 1);
    fplot(a_post*sin(t-pi/2), diff(k_elip)^2, elip_bounds, 'LineWidth', 1);
    fplot(a_post_fourier*sin(t-pi), diff(k_fourier)^2, fourier_bounds, 'LineWidth', 1);
    fplot(t, diff(k_forbes)^2, 'LineWidth', 1);
    set(gca,'TickDir','out'); ax=gca; ax.FontSize=16;
    legend("Chien", "Ellipse", "Fourier", "Forbes", 'FontSize', 20); title("Square of the First Derivative of the Posterior Curvature for the 20-Year Old Lens", 'FontSize', 24); xlabel("Optic Zone Radius (mm)", 'FontSize', 20); ylabel("First Derivative of Curvature Squared (1/mm^4)", 'FontSize', 20); ylim([0, 0.8]); xlim([-3,3]);
    figure
    
%     curve_title = ["Raw X", "Raw Y", "Chien X", "Chien Y", "Ellipse X", "Ellipse Y", "Fourier X", "Fourier Y", "Forbes X", "Forbes Y"];
%     curve_data = [X_data'; Y_data'; X_chien,NaN(1,1413); Y_chien,NaN(1,1413); X_elip,NaN(1,1431); Y_elip,NaN(1,1431); X_fourier,NaN(1,1422); Y_fourier,NaN(1,1422); X_forbes,NaN(1,1524); Y_forbes,NaN(1,1524)];
%     k_title = ["Curvature Chien X", "Curvature Chien Y", "Curvature Ellipse X", "Curvature Ellipse Y", "Curvature Fourier X", "Curvature Fourier Y", "Curvature Forbes X", "Curvature Forbes Y"];
%     k_data = [K_chien_x,NaN(1,398); K_chien_y,NaN(1,398); K_elip_x,NaN(1,510); K_elip_y,NaN(1,510); K_fourier_x; K_fourier_y; K_forbes_x,NaN(1,493); K_forbes_y,NaN(1,493)];
%     
%     writematrix(curve_title, "post_graph_data.xls", "Sheet", "Posterior Raw Data", "Range", 'A1')
%     writematrix(curve_data', "post_graph_data.xls", "Sheet", "Posterior Raw Data", "Range", 'A2')
%     
%     writematrix(k_title, "post_graph_data.xls", "Sheet", "Posterior Curvature Data", "Range", 'A1')
%     writematrix(k_data', "post_graph_data.xls", "Sheet", "Posterior Curvature Data", "Range", 'A2')
    
end


