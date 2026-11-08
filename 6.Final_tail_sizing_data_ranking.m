%% DEFINITIVE PARAMETRIC TAIL RANKING PIPELINE (FLOW5 VLM2)
% Department of Mechanical Engineering, University of Moratuwa
clear; clc; close all;

zip_filename = 'Tail_analysis_201.zip';
desktop_path = fullfile(getenv('USERPROFILE'), 'Desktop', zip_filename);
unzip_target = 'Tail_analysis_201';

% 1. Extract zip archive if not already unpacked
if ~exist(unzip_target, 'dir')
    if exist(zip_filename, 'file')
        fprintf('Extracting %s...\n', zip_filename);
        unzip(zip_filename, unzip_target);
    elseif exist(desktop_path, 'file')
        fprintf('Extracting %s from Desktop...\n', zip_filename);
        unzip(desktop_path, unzip_target);
    else
        error('Cannot find Tail_analysis_201.zip or extracted folder.');
    end
else
    fprintf('Using extracted folder "%s".\n', unzip_target);
end

% Glider Reference Geometry (from flow5 Plane Explorer)
Sw    = 0.024; % Main wing reference area (m^2)
c_mac = 0.060; % Mean Aerodynamic Chord (m)

% 2. Locate all polar text files
files = dir(fullfile(unzip_target, '**', '*.txt'));
if isempty(files)
    files = dir(fullfile(unzip_target, '**', '*.*'));
    valid_idx = ~[files.isdir] & ~endsWith({files.name}, {'.m', '.csv', '.zip', '.png'});
    files = files(valid_idx);
end

num_files = length(files);
fprintf('Found %d tail polar files. Parsing aerodynamic sweeps...\n', num_files);

tail_results = struct();
valid_count = 0;

% 3. Robust In-Memory Parser
for k = 1:num_files
    fp = fullfile(files(k).folder, files(k).name);
    fid = fopen(fp, 'r');
    if fid == -1, continue; end
    raw_text = fread(fid, '*char')';
    fclose(fid);
    
    lines = splitlines(raw_text);
    num_lines = length(lines);
    if num_lines < 6, continue; end
    
    alpha = []; CL = []; CD = []; Cm = []; CL_CD = []; CL15_CD = [];
    
    for L = 1:num_lines
        line_str = strtrim(lines{L});
        if isempty(line_str), continue; end
        
        % In flow5 exports, Line 1 of each operating point begins with "0"
        if line_str(1) == '0'
            vals = sscanf(line_str, '%f');
            
            % Line 1 contains: Ctrl, Alpha, Beta, Phi, CL, CD, ..., Cm (Col 10)
            if length(vals) >= 10 && abs(vals(1)) < 1e-4
                a_val  = vals(2);
                cl_val = vals(5);
                cd_val = vals(6);
                cm_val = vals(10);
                
                % Filter unphysical drag zeros and compute ratios directly
                if cd_val > 0.005 && cl_val > 0.01
                    alpha(end+1, 1)   = a_val;
                    CL(end+1, 1)      = cl_val;
                    CD(end+1, 1)      = cd_val;
                    Cm(end+1, 1)      = cm_val;
                    CL_CD(end+1, 1)   = cl_val / cd_val;
                    CL15_CD(end+1, 1) = (cl_val^1.5) / cd_val;
                end
            end
        end
    end
    
    if length(alpha) < 3, continue; end
    
    valid_count = valid_count + 1;
    [~, cfg_name, ~] = fileparts(files(k).name);
    
    % Geometry Extraction (Handles decimal underscores: 4_5 -> 4.5, 3_5 -> 3.5)
    tok = regexp(cfg_name, '^(\d+)_(\d+)_(\d+)_(\d+)', 'tokens');
    if ~isempty(tok)
        St_val = str2double(tok{1}{1});
        AR_val = str2double([tok{1}{2} '.' tok{1}{3}]);
        lt_val = str2double(tok{1}{4});
    else
        tok2 = regexp(cfg_name, '^(\d+)_([\d\.]+)_(\d+)', 'tokens');
        if ~isempty(tok2)
            St_val = str2double(tok2{1}{1});
            AR_val = str2double(tok2{1}{2});
            lt_val = str2double(tok2{1}{3});
        else
            St_val = NaN; AR_val = NaN; lt_val = NaN;
        end
    end
    
    % True Vh calculation
    if ~isnan(St_val) && ~isnan(lt_val)
        St_m2  = St_val * 1e-6;
        lt_m   = lt_val * 1e-3;
        Vh_val = (St_m2 * lt_m) / (Sw * c_mac);
    else
        Vh_val = NaN;
    end
    
    [max_ld, idx_ld] = max(CL_CD);
    max_pf           = max(CL15_CD);
    cd_min           = min(CD);
    
    p = polyfit(alpha, Cm, 1);
    cm_alpha = p(1);
    
    tail_results(valid_count).ConfigName   = cfg_name;
    tail_results(valid_count).St           = St_val;
    tail_results(valid_count).ARt          = AR_val;
    tail_results(valid_count).lt           = lt_val;
    tail_results(valid_count).Vh           = Vh_val;
    tail_results(valid_count).Max_LD       = max_ld;
    tail_results(valid_count).Alpha_Opt    = alpha(idx_ld);
    tail_results(valid_count).CL_at_MaxLD  = CL(idx_ld);
    tail_results(valid_count).CD_at_MaxLD  = CD(idx_ld);
    tail_results(valid_count).CD_min       = cd_min;
    tail_results(valid_count).Max_PF       = max_pf;
    tail_results(valid_count).Cm_alpha     = cm_alpha;
    tail_results(valid_count).alpha        = alpha;
    tail_results(valid_count).CL_CD        = CL_CD;
    tail_results(valid_count).Cm           = Cm;
end

if valid_count == 0
    error('Could not extract valid operating polars. Verify polar data format.');
end

%% 4. Multi-Attribute Scoring & Flight Stability Filtering
all_ld   = [tail_results.Max_LD];
all_pf   = [tail_results.Max_PF];
all_cm_a = [tail_results.Cm_alpha];

% Feature Normalization (0.0 to 1.0)
norm_ld = (all_ld - min(all_ld)) ./ (max(all_ld) - min(all_ld) + eps);
norm_pf = (all_pf - min(all_pf)) ./ (max(all_pf) - min(all_pf) + eps);

stiff_val = -all_cm_a;
norm_cm   = (stiff_val - min(stiff_val)) ./ (max(stiff_val) - min(stiff_val) + eps);

for i = 1:valid_count
    % 50% Glide Efficiency + 35% Flight Time + 15% Pitch Stiffness
    base_score = (0.50 * norm_ld(i) + 0.35 * norm_pf(i) + 0.15 * norm_cm(i)) * 100;
    
    vh_curr = tail_results(i).Vh;
    % Practical free-flight stability envelope: 0.40 <= Vh <= 0.60
    if ~isnan(vh_curr) && (vh_curr < 0.40 || vh_curr > 0.60)
        tail_results(i).Status = 'Outside Vh [0.4-0.6]';
        tail_results(i).Score  = base_score * 0.65;
    else
        tail_results(i).Status = 'Qualified Stable';
        tail_results(i).Score  = base_score;
    end
end

[~, sort_idx] = sort([tail_results.Score], 'descend');
sorted = tail_results(sort_idx);

%% 5. Print Output Table
fprintf('\n========================================================================================================================================\n');
fprintf('%-5s %-32s %-7s %-10s %-10s %-9s %-18s %-10s %-7s %-18s\n', ...
    'Rank', 'Configuration', 'Vh', '3D Max L/D', 'Opt Alpha', 'CD_min', 'Cm_alpha (deg^-1)', 'CL^1.5/CD', 'Score', 'Status');
fprintf('========================================================================================================================================\n');
for r = 1:min(20, valid_count)
    fprintf('%-5d %-32s %-7.3f %-10.2f %-10.1f %-9.4f %-18.4f %-10.2f %-7.1f %-18s\n', ...
        r, sorted(r).ConfigName(1:min(32, length(sorted(r).ConfigName))), ...
        sorted(r).Vh, sorted(r).Max_LD, sorted(r).Alpha_Opt, sorted(r).CD_min, ...
        sorted(r).Cm_alpha, sorted(r).Max_PF, sorted(r).Score, sorted(r).Status);
end
fprintf('========================================================================================================================================\n\n');

w = sorted(1);
fprintf('>>> WINNING TAIL CONFIGURATION VERDICT <<<\n');
fprintf('  Configuration Name    : %s\n', w.ConfigName);
fprintf('  Tail Volume (Vh)      : %.3f\n', w.Vh);
fprintf('  Horizontal Area (St)  : %.1f mm^2\n', w.St);
fprintf('  Tail Aspect Ratio     : %.2f\n', w.ARt);
fprintf('  Tail Moment Arm (lt)  : %.1f mm\n', w.lt);
fprintf('  Peak Finite Glide L/D : %.2f (at alpha = %.1f deg)\n', w.Max_LD, w.Alpha_Opt);
fprintf('  Parasite Drag Floor   : CD_min = %.4f\n', w.CD_min);
fprintf('  Pitch Stiffness       : Cm_alpha = %.4f deg^-1 (%.3f rad^-1)\n', w.Cm_alpha, w.Cm_alpha * (180/pi));
fprintf('  Composite Score       : %.1f / 100\n\n', w.Score);

%% 6. Export Results to CSV
csv_name = 'TAIL_ANALYSIS_RANKED.csv';
fid = fopen(csv_name, 'w');
if fid ~= -1
    fprintf(fid, 'Rank,ConfigName,St_mm2,ARt,lt_mm,Vh,Max_LD,Alpha_Opt,CL_at_MaxLD,CD_at_MaxLD,CD_min,Max_CL15_CD,Cm_alpha_deg,Cm_alpha_rad,Score,Status\n');
    for r = 1:valid_count
        it = sorted(r);
        fprintf(fid, '%d,%s,%.1f,%.2f,%.1f,%.3f,%.2f,%.2f,%.3f,%.4f,%.4f,%.2f,%.4f,%.3f,%.1f,%s\n', ...
            r, it.ConfigName, it.St, it.ARt, it.lt, it.Vh, it.Max_LD, it.Alpha_Opt, ...
            it.CL_at_MaxLD, it.CD_at_MaxLD, it.CD_min, it.Max_PF, it.Cm_alpha, ...
            it.Cm_alpha * (180/pi), it.Score, it.Status);
    end
    fclose(fid);
    fprintf('[OK] Full ranked matrix exported to "%s"\n', csv_name);
end

%% 7. Generate Comparison Plots
num_plot = min(5, valid_count);
if num_plot > 0
    fig = figure('Name', 'Optimal Tail Configurations Comparison', 'Color', 'w', 'Position', [100, 100, 1100, 480]);
    plot_colors = {'b', 'r', [0 0.6 0], 'm', 'k'};

    subplot(1, 2, 1); hold on; grid on; box on;
    for i = 1:num_plot
        c_use = plot_colors{mod(i-1, length(plot_colors)) + 1};
        plot(sorted(i).alpha, sorted(i).CL_CD, 'LineWidth', 1.8, 'Color', c_use, ...
            'DisplayName', sprintf('Rank %d: %s (Vh=%.2f)', i, sorted(i).ConfigName, sorted(i).Vh));
    end
    xlabel('Angle of Attack \alpha (deg)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Glide Ratio (C_L / C_D)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Finite Glide Efficiency Comparison', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'none');

    subplot(1, 2, 2); hold on; grid on; box on;
    for i = 1:num_plot
        c_use = plot_colors{mod(i-1, length(plot_colors)) + 1};
        plot(sorted(i).alpha, sorted(i).Cm, 'LineWidth', 1.8, 'Color', c_use, ...
            'DisplayName', sprintf('%s (Cm_\\alpha=%.3f)', sorted(i).ConfigName, sorted(i).Cm_alpha));
    end
    xlabel('Angle of Attack \alpha (deg)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Pitching Moment Coefficient (C_m)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Pitch Restoring Stability Slope (C_m vs \alpha)', 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location', 'northeast', 'FontSize', 8, 'Interpreter', 'none');

    saveas(fig, 'Figure_Tail_Optimization_Curves.png');
    fprintf('[OK] Performance plots saved as "Figure_Tail_Optimization_Curves.png"\n');
end