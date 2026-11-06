%% DEFINITIVE SYMMETRIC TAIL AIRFOIL RANKING & COMPARISON SCRIPT
clear; clc; close all;

zip_name = 'Selected_tail_aerofoil_polars.zip';
data_dir = 'Selected_tail_aerofoil_polars';

% Auto-unzip if folder is not already extracted
if ~exist(data_dir, 'dir')
    if exist(zip_name, 'file')
        fprintf('Extracting %s...\n', zip_name);
        unzip(zip_name, data_dir);
    else
        error('Neither folder "%s" nor zip file "%s" was found in current directory.', data_dir, zip_name);
    end
end

% Recursively gather all polar text files
files = dir(fullfile(data_dir, '**', '*.*'));
valid_exts = {'.txt', '.plr', '.csv', '.dat', ''};

file_list = {};
for i = 1:length(files)
    [~, ~, ext] = fileparts(files(i).name);
    if ~files(i).isdir && ismember(lower(ext), valid_exts) && ...
       ~contains(files(i).name, 'RANKING') && ~endsWith(files(i).name, '.m')
        file_list{end+1} = fullfile(files(i).folder, files(i).name); %#ok<SAGROW>
    end
end

if isempty(file_list)
    error('No polar files found inside "%s".', data_dir);
end

tail_map = containers.Map();

%% 1. Parse Each Polar File
for k = 1:length(file_list)
    fp = file_list{k};
    fid = fopen(fp, 'r');
    if fid == -1, continue; end
    
    foil_name = '';
    re_str    = 'N/A';
    header_lines = 0;
    
    % Read headers to extract Foil Name, Reynolds Number, and header depth
    while ~feof(fid)
        tline = fgetl(fid);
        header_lines = header_lines + 1;
        if ~ischar(tline), break; end
        
        if contains(tline, 'Calculated polar for:')
            raw_name = strtrim(strrep(tline, 'Calculated polar for:', ''));
            % Clean trailing thickness annotations like '(15.01%)'
            foil_name = strtrim(regexprep(raw_name, '\(.*?\)', ''));
        elseif contains(tline, 'Re =')
            tok = regexp(tline, 'Re\s*=\s*([0-9\.\s\w\+\-]+?)(?:\s+Ncrit|\s*$)', 'tokens');
            if ~isempty(tok)
                re_raw = strtrim(tok{1}{1});
                % Convert scientific format '0.002 e 6' -> '2,000'
                num_tok = regexp(re_raw, '([0-9\.]+)\s*[eE]\s*([0-9\+]+)', 'tokens');
                if ~isempty(num_tok)
                    val = str2double(num_tok{1}{1}) * (10^str2double(num_tok{1}{2}));
                    re_str = sprintf('%d', round(val));
                else
                    re_str = re_raw;
                end
            end
        end
        
        % Detect end of header (the line with dashes below column names)
        if contains(tline, '-----') && header_lines >= 8
            break;
        end
    end
    fclose(fid);
    
    if isempty(foil_name)
        [~, base_name, ~] = fileparts(fp);
        foil_name = regexprep(base_name, '_T\d+.*$', '');
    end
    
    % Read numerical polar columns
    try
        mat = readmatrix(fp, 'HeaderLines', header_lines);
    catch
        continue;
    end
    
    if isempty(mat) || size(mat, 1) < 2 || size(mat, 2) < 3
        continue;
    end
    
    alpha = mat(:, 1);
    cl    = mat(:, 2);
    cd    = mat(:, 3);
    
    valid = ~isnan(alpha) & ~isnan(cl) & ~isnan(cd);
    alpha = alpha(valid);
    cl    = cl(valid);
    cd    = cd(valid);
    
    if length(alpha) < 2, continue; end
    
    % Zero-lift parasite drag at alpha closest to 0 deg
    [~, idx_0] = min(abs(alpha - 0.0));
    cd0_val = cd(idx_0);
    
    % Lift curve slope near neutral pitch [-2.0 deg, +2.0 deg]
    fit_idx = (alpha >= -2.0) & (alpha <= 2.0);
    if sum(fit_idx) >= 2
        p = polyfit(alpha(fit_idx), cl(fit_idx), 1);
        dcl_da = p(1);
    else
        dcl_da = (cl(end) - cl(1)) / (alpha(end) - alpha(1));
    end
    
    % Store unique profile data (retain lower Cd0 if duplicate runs exist)
    if ~isKey(tail_map, foil_name) || cd0_val < tail_map(foil_name).Cd0
        entry.Name    = foil_name;
        entry.Re      = re_str;
        entry.Cd0     = round(cd0_val, 4);
        entry.dCl_da  = round(dcl_da, 4);
        entry.dCl_rad = round(dcl_da * (180 / pi), 3);
        entry.alpha   = alpha;
        entry.cl      = cl;
        entry.cd      = cd;
        tail_map(foil_name) = entry;
    end
end

keys = tail_map.keys();
if isempty(keys)
    error('Could not extract valid data points from polar files.');
end

N = length(keys);
records = cell(N, 1);
for i = 1:N
    records{i} = tail_map(keys{i});
end

%% 2. Multi-Attribute Scoring (50% Slope Authority + 50% Low Drag Floor)
all_slopes = cellfun(@(x) x.dCl_da, records);
all_cd0    = cellfun(@(x) x.Cd0, records);

norm_slope = (all_slopes - min(all_slopes)) ./ (max(all_slopes) - min(all_slopes) + eps);
norm_drag  = 1.0 - ((all_cd0 - min(all_cd0)) ./ (max(all_cd0) - min(all_cd0) + eps));

for i = 1:N
    records{i}.Score = (0.50 * norm_slope(i) + 0.50 * norm_drag(i)) * 100;
end

[~, sort_idx] = sort(cellfun(@(x) x.Score, records), 'descend');
sorted_tails = records(sort_idx);

%% 3. Print Complete Ranked Leaderboard to Console
fprintf('\n==================================================================================================================\n');
fprintf('Rank  Airfoil Name             Re          Score   Cd0 (alpha=0)   dCl/da (deg^-1)  dCl/da (rad^-1)  Deadband Rating\n');
fprintf('==================================================================================================================\n');
for i = 1:N
    t = sorted_tails{i};
    if t.dCl_da >= 0.085
        rating = 'Excellent (Sharp)';
    elseif t.dCl_da >= 0.070
        rating = 'Good';
    elseif t.dCl_da >= 0.040
        rating = 'Moderate';
    else
        rating = 'Poor (Severe Deadband)';
    end
    fprintf('%-5d %-24s %-11s %-7.1f %-15.4f %-16.4f %-16.3f %-15s\n', ...
        i, t.Name(1:min(23, length(t.Name))), t.Re, t.Score, t.Cd0, t.dCl_da, t.dCl_rad, rating);
end
fprintf('==================================================================================================================\n\n');

fprintf('>>> WINNING TAIL AIRFOIL VERDICT <<<\n');
fprintf('  Selected Profile     : %s\n', sorted_tails{1}.Name);
fprintf('  Parasite Drag at 0°  : Cd0 = %.4f\n', sorted_tails{1}.Cd0);
fprintf('  Lift Slope near 0°   : dCl/da = %.4f deg^-1 (%.3f rad^-1)\n', sorted_tails{1}.dCl_da, sorted_tails{1}.dCl_rad);
fprintf('  Overall Tail Score   : %.1f / 100\n\n', sorted_tails{1}.Score);

%% 4. Generate Comparative Aerodynamic Plots
fig = figure('Name', 'Symmetric Tail Airfoil Evaluation', 'Color', 'w', 'Position', [80, 80, 1100, 500]);
clf(fig);
palette = lines(N);

% Subplot 1: Lift Slope & Control Deadband Check (Cl vs Alpha)
subplot(1, 2, 1);
hold on; grid on; box on;
for i = 1:N
    t = sorted_tails{i};
    [a_sort, s_idx] = sort(t.alpha);
    plot(a_sort, t.cl(s_idx), '-o', 'Color', palette(i, :), 'LineWidth', 2.0, ...
        'MarkerSize', 4, 'DisplayName', sprintf('%s (Score: %.1f)', t.Name, t.Score));
end
xlim([-4, 4]);
xlabel('Angle of Attack \alpha (deg)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Section Lift Coefficient (C_l)', 'FontSize', 11, 'FontWeight', 'bold');
title('Lift Curve & Deadband Sensitivity (\alpha = 0°)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 8, 'Interpreter', 'none');

% Subplot 2: Parasite Drag Floor (Cd vs Alpha)
subplot(1, 2, 2);
hold on; grid on; box on;
for i = 1:N
    t = sorted_tails{i};
    [a_sort, s_idx] = sort(t.alpha);
    plot(a_sort, t.cd(s_idx), '-s', 'Color', palette(i, :), 'LineWidth', 2.0, ...
        'MarkerSize', 4, 'DisplayName', t.Name);
end
xlim([-4, 4]);
xlabel('Angle of Attack \alpha (deg)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Section Drag Coefficient (C_d)', 'FontSize', 11, 'FontWeight', 'bold');
title('Tail Parasite Drag Floor (C_d vs \alpha)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'north', 'FontSize', 8, 'Interpreter', 'none');

drawnow;
saveas(fig, 'Figure_5_Tail_Airfoil_Comparison.png');
fprintf('[OK] Figure saved as "Figure_5_Tail_Airfoil_Comparison.png"\n');

%% 5. Export Ranked CSV Table
csv_out = 'TAIL_AIRFOILS_RANKED.csv';
fid = fopen(csv_out, 'w');
if fid ~= -1
    fprintf(fid, 'Rank,Airfoil,Re,Score,Cd0,dCl_da_deg,dCl_da_rad,Rating\n');
    for i = 1:N
        t = sorted_tails{i};
        if t.dCl_da >= 0.085
            r_str = 'Excellent (Sharp)';
        elseif t.dCl_da >= 0.070
            r_str = 'Good';
        elseif t.dCl_da >= 0.040
            r_str = 'Moderate';
        else
            r_str = 'Poor (Severe Deadband)';
        end
        fprintf(fid, '%d,"%s","%s",%.2f,%.4f,%.4f,%.3f,"%s"\n', ...
            i, t.Name, t.Re, t.Score, t.Cd0, t.dCl_da, t.dCl_rad, r_str);
    end
    fclose(fid);
    fprintf('[OK] Table exported to: %s\n', csv_out);
end