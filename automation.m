clear; clc; close all;
% Main script. Run this to replicate results %

%% Definitions
% Data files to do analysis on
data_file = {'Human_20 yo RIEB15-1632_OD_data.xls', 'Human_21 yo RIEB14-1748_OD_data.xls', 'Human_22 yo RIEB13-1595_OD_data.xls', 'Human_24 yo RIEB13-0253_OS_data', 'Human_25 yo RIEB15-1976_OD_data.xls', 'Human_26 yo RIEB14-1243_OD_data.xls', 'Human_27 yo RIEB16-0368_OD_data.xls', 'Human_28 yo RIEB13-1936_OD_data.xls', 'Human_29 yo RIEB13-0768T1_OD_data.xls', 'Human_30 yo RIEB13-0161_OD_data.xls'};

% For file naming purposes
age_vec = [20:22,24:30];
chien_curve = "5";
forbes_conic = "elip-regWeight";

%% Generate data template
chien_label = strcat('Chien (Curve ', chien_curve, ')');
forbes_label = strcat('Forbes (Conic=', forbes_conic, ')');

label_cell = cell(13,11);
label_cell(1,:) = {'Anterior Data', chien_label, forbes_label,'Fourier', 'ellipse','', 'Posterior Data', chien_label, forbes_label,'Fourier', 'ellipse'};
label_cell(2:end,1) = {'Fit', 'Arc Length', 'Vertex Curvature', 'Bending Energy', 'Waviness', 'Variance of Curvature', 'ARoC 1', 'ARoC 2', 'ARoC 3', 'ARoC 4', 'ARoC 5', 'ARoC 6'};
label_cell(2:end,7)={'Fit', 'Arc Length', 'Vertex Curvature', 'Bending Energy', 'Waviness', 'Variance of Curvature', 'PRoC 1', 'PRoC 2', 'PRoC 3', 'PRoC 4', 'PRoC 5', 'PRoC 6'};

%% Get results
% For each file
for i = 1:length(data_file)
    file = data_file(i)
    age = age_vec(i);
    
    full_cell=label_cell;
        
    % Get anterior results
    M=0;
    ant_results = anterior(file, M);
    close all;
    
    % Get posterior results
    post_results = posterior(file, M);
    close all;
    
    % Store results
    full_cell(2:end, 2:5) = num2cell(ant_results);
    full_cell(2:end, 8:end) = num2cell(post_results);
    
    % Save results in Excel sheet
    writecell(full_cell, "results.xls", 'sheet', strcat('age', num2str(age)))
end