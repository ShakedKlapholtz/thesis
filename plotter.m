%% PARAMETERS
N = 256;
cascade_block = 73;
I_4 =6.47;
p = 0.01; % Define the probability for the cascade

file_path_csv = "\\wsl.localhost\Ubuntu-22.04\home\shaked\shake-on-it\cascade_reconcile_p_0.10_256_runs.csv";
field_name    = 'reply_parity_bits';

file_path_fig = "C:\Users\ASUS\Documents\MATLAB\Key_Length_256_SCL_L=32_-0.8_100_trials_BSC.fig";

%% ===============================
%  LOAD CASCADE DATA (CSV)
% ===============================
T = readtable(file_path_csv);

if ~ismember(field_name, T.Properties.VariableNames)
    error('Field %s not found in CSV file', field_name);
end

% Exposure = N - key length
Exposure_cascade = T.(field_name);
Exposure_cascade_norm = Exposure_cascade / N;

%% ===============================
%  LOAD POLAR HISTOGRAM (FIG)
% ===============================
fig = openfig(file_path_fig, 'invisible');
ax = gca(fig);

h_existing = findobj(ax, 'Type', 'histogram');
if isempty(h_existing)
    close(fig);
    error('No histogram found in FIG file');
end

existingCounts = h_existing.BinCounts;
existingEdges  = N - h_existing.BinEdges;

close(fig);

% Normalize Polar histogram
existingProb = existingCounts / sum(existingCounts);
binCenters   = (existingEdges(1:end-1) + existingEdges(2:end)) / 2;

% Convert to relative exposure
binCenters_norm = binCenters / N;

%% ===============================
%  HISTOGRAM SETTINGS
% ===============================
edges_norm = linspace(0,1,N+1);

%% ===============================
%  MEANS
% ===============================
mean_polar_norm   = sum(binCenters_norm .* existingProb);
mean_cascade_norm = mean(Exposure_cascade_norm);
mean_cascade_theory = (p*I_4 + 1/cascade_block);

%% ===============================
%  PLOT
% ===============================
figure;
hold on;

% Polar (from FIG)
bar(binCenters_norm, existingProb, ...
    'FaceColor', [0.8 0.8 1], ...
    'EdgeColor', 'b', ...
    'DisplayName', 'Polar');

% Cascade (from CSV)
histogram(Exposure_cascade_norm, edges_norm, ...
    'Normalization', 'probability', ...
    'FaceColor', [1 0.8 0.8], ...
    'FaceAlpha', 0.6, ...
    'EdgeColor', 'r', ...
    'DisplayName', 'Cascade');

% Mean lines
xline(mean_polar_norm, '--b', ...
    sprintf('Mean Polar = %.3f', mean_polar_norm), ...
    'LineWidth', 2, 'DisplayName', 'Mean Polar');

xline(mean_cascade_norm, '--r', ...
    sprintf('Mean Cascade = %.3f', mean_cascade_norm), ...
    'LineWidth', 2, 'DisplayName', 'Mean Cascade');

xline(mean_cascade_theory, '--g', ...
    sprintf('Mean Cascade Theory = %.3f', mean_cascade_theory), ...
    'LineWidth', 2, 'DisplayName', 'Mean Cascade in theory');

xlabel('Relative Exposure');
ylabel('Probability');
title(sprintf('Relative Exposure Distribution (N=%d, p=%.2f)', N, p));
legend show;
grid on;
hold off;

%% ===============================
%  DISPLAY MEANS
% ===============================
fprintf('Mean relative exposure (Polar)   = %.4f\n', mean_polar_norm);
fprintf('Mean relative exposure (Cascade) = %.4f\n', mean_cascade_norm);
fprintf('Mean relative exposure in theory (Cascade) = %.4f\n', mean_cascade_theory);

