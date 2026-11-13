% =========================================================================
% FINAL FULL-GLIDER 3D AERODYNAMIC, TRIM & STABILITY VERIFICATION
% =========================================================================
clear; clc; close all;

%% 1. Load Polar File Directly from Current Folder
targetFile = 'Final_Glider_Design_T2-VLM2-0_0kg-x0_0m.txt';
if ~exist(targetFile, 'file')
    txtFiles = dir('*.txt');
    if isempty(txtFiles)
        error('No .txt file found in the current MATLAB folder.');
    end
    targetFile = txtFiles(1).name;
end

m   = 0.035;        % Total Design Mass (kg)
g   = 9.81;         % Gravity (m/s^2)
W   = m * g;        % Design Weight (0.34335 N)
rho = 1.225;        % Air Density (kg/m^3)
S_w = 0.02353;      % Main Wing Area (m^2)

%% 2. Parse Primary Data Rows (Single-Line Alignment — Zero Array Mismatch)
fprintf('Reading Polar File: %s\n', targetFile);
lines = splitlines(fileread(targetFile));

alpha_raw = [];
CL_raw    = [];
CD_raw    = [];
Cm_raw    = [];

for k = 1:length(lines)
    tline = strtrim(lines{k});
    if isempty(tline), continue; end
    
    % Primary data rows in flow5 VLM2 have >= 14 numeric columns starting with Ctrl = 0
    v = sscanf(tline, '%f');
    if length(v) >= 14 && abs(v(1)) < 1e-5 && abs(v(3)) < 1e-5 && abs(v(4)) < 1e-5
        alpha_raw(end+1, 1) = v(2);  %#ok<AGROW> % Column 2: Alpha (deg)
        CL_raw(end+1, 1)    = v(5);  %#ok<AGROW> % Column 5: 3D Lift (CL)
        CD_raw(end+1, 1)    = v(6);  %#ok<AGROW> % Column 6: 3D Total Drag (CD)
        Cm_raw(end+1, 1)    = v(10); %#ok<AGROW> % Column 10: Pitching Moment (Cm)
    end
end

if length(alpha_raw) < 4
    error('Failed to extract polar rows. Verify file format.');
end

% Sort sequentially by Angle of Attack (all arrays guaranteed identical length)
[alpha_raw, idx_sort] = sort(alpha_raw);
CL_raw = CL_raw(idx_sort);
CD_raw = CD_raw(idx_sort);
Cm_raw = Cm_raw(idx_sort);
LD_raw = CL_raw ./ CD_raw;

%% 3. High-Resolution Spline Interpolation
alpha_fine = linspace(min(alpha_raw), max(alpha_raw), 2000)';
CL_fine    = spline(alpha_raw, CL_raw, alpha_fine);
CD_fine    = spline(alpha_raw, CD_raw, alpha_fine);
Cm_fine    = spline(alpha_raw, Cm_raw, alpha_fine);
LD_fine    = CL_fine ./ CD_fine;

[max_LD, idx_LD] = max(LD_fine);
alpha_max_LD     = alpha_fine(idx_LD);
CL_max_LD        = CL_fine(idx_LD);
CD_max_LD        = CD_fine(idx_LD);

%% 4. Longitudinal Trim Equilibrium (Cm = 0)
sign_changes = find(diff(sign(Cm_fine)) ~= 0);
trim_alpha   = NaN;
CL_trim      = NaN;
CD_trim      = NaN;
LD_trim      = NaN;
V_trim       = NaN;

if ~isempty(sign_changes)
    idx = sign_changes(1);
    x1 = alpha_fine(idx); x2 = alpha_fine(idx+1);
    y1 = Cm_fine(idx);    y2 = Cm_fine(idx+1);
    trim_alpha = x1 - y1 * (x2 - x1) / (y2 - y1);
    
    CL_trim = interp1(alpha_fine, CL_fine, trim_alpha, 'spline');
    CD_trim = interp1(alpha_fine, CD_fine, trim_alpha, 'spline');
    LD_trim = CL_trim / CD_trim;
    if CL_trim > 0
        V_trim = sqrt((2 * W) / (rho * S_w * CL_trim));
    end
end

%% 5. Static Pitch Stability Slope (Cm_alpha)
lin_mask     = alpha_raw <= 7.5;
poly_fit     = polyfit(alpha_raw(lin_mask), Cm_raw(lin_mask), 1);
Cm_alpha_deg = poly_fit(1);
Cm_alpha_rad = Cm_alpha_deg * (180 / pi);

if Cm_alpha_deg < 0
    stability_verdict = 'STATICALLY STABLE (Inherent Nose-Down Restoring Torque)';
else
    stability_verdict = 'STATICALLY UNSTABLE (Divergent Pitch Moment)';
end

%% 6. Print Complete Verification Report
fprintf('\n=====================================================================\n');
fprintf('        FINAL FULL-GLIDER 3D AERODYNAMIC & STABILITY REPORT          \n');
fprintf('=====================================================================\n');
fprintf(' Stability Verdict          : %s\n', stability_verdict);
fprintf(' Pitch Stiffness (Cm_alpha) : %.4f deg^-1  (%.4f rad^-1)\n', Cm_alpha_deg, Cm_alpha_rad);
fprintf('---------------------------------------------------------------------\n');
fprintf(' Equilibrium Trim Angle     : %.3f deg  (at Cm = 0)\n', trim_alpha);
fprintf(' Trim Lift Coefficient (CL) : %.3f\n', CL_trim);
fprintf(' Trim Drag Coefficient (CD) : %.4f\n', CD_trim);
fprintf(' Trim Glide Ratio (L/D)     : %.2f\n', LD_trim);
fprintf(' Equilibrium Glide Speed    : %.2f m/s\n', V_trim);
fprintf('---------------------------------------------------------------------\n');
fprintf(' Peak 3D Glide Ratio (L/D)  : %.2f  (at alpha = %.3f deg)\n', max_LD, alpha_max_LD);
fprintf(' Peak L/D Lift & Drag       : CL = %.3f, CD = %.4f\n', CL_max_LD, CD_max_LD);
fprintf('=====================================================================\n\n');

%% 7. Generate & Export Verification Graphs
fig = figure('Name', 'Full Glider Verification Report', 'Color', 'w', 'Position', [80, 100, 1150, 480]);

subplot(1, 2, 1); hold on; grid on; box on;
plot(alpha_raw, Cm_raw, 'ko', 'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerSize', 6, ...
    'DisplayName', 'flow5 VLM2 Data');
plot(alpha_fine, Cm_fine, 'b-', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('C_m Curve (C_{m\\alpha} = %.4f deg^{-1})', Cm_alpha_deg));
yline(0, 'k-', 'LineWidth', 1.1, 'HandleVisibility', 'off');

if ~isnan(trim_alpha)
    plot(trim_alpha, 0, 'r*', 'MarkerSize', 12, 'LineWidth', 1.6, ...
        'DisplayName', sprintf('Trim Angle: %.3f°', trim_alpha));
    xline(trim_alpha, 'r--', 'HandleVisibility', 'off');
end

xlabel('Angle of Attack, \alpha (deg)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Pitching Moment Coefficient (C_m)', 'FontSize', 11, 'FontWeight', 'bold');
title('Pitch Stability (C_m vs \alpha)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 9);

subplot(1, 2, 2); hold on; grid on; box on;
plot(alpha_raw, LD_raw, 'ko', 'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerSize', 6, ...
    'DisplayName', 'flow5 VLM2 Data');
plot(alpha_fine, LD_fine, 'g-', 'LineWidth', 2.0, ...
    'DisplayName', 'Efficiency Curve');
plot(alpha_max_LD, max_LD, 'p', 'MarkerSize', 15, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'r', ...
    'LineWidth', 1.2, 'DisplayName', sprintf('Max L/D: %.2f (at %.2f°)', max_LD, alpha_max_LD));

if ~isnan(trim_alpha)
    plot(trim_alpha, LD_trim, 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'c', ...
        'DisplayName', sprintf('Trim L/D: %.2f (at %.2f°)', LD_trim, trim_alpha));
end

xlabel('Angle of Attack, \alpha (deg)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Efficiency (C_L / C_D)', 'FontSize', 11, 'FontWeight', 'bold');
title('Aerodynamic Efficiency (C_L/C_D vs \alpha)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'south', 'FontSize', 8.5);

saveas(fig, 'Figure_Full_Glider_Verification.png');
fprintf('[OK] Plot saved as "Figure_Full_Glider_Verification.png"\n');