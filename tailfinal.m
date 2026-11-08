%% PARAMETRIC TAIL RANKING WITH DESIGN-FOR-MANUFACTURING (DFM)
% Department of Mechanical Engineering, University of Moratuwa
clear; clc; close all;

zip_filename = 'Tail_analysis_201.zip';
desktop_path = fullfile(getenv('USERPROFILE'), 'Desktop', zip_filename);
unzip_target = 'Tail_analysis_201';

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

% Glider Reference Geometry
Sw    = 0.024; % Main wing reference area (m^2)
c_mac = 0.060; % Mean Aerodynamic Chord (m)

files = dir(fullfile(unzip_target, '**', '*.txt'));
if isempty(files)
    files = dir(fullfile(unzip_target, '**', '*.*'));
    valid_idx = ~[files.isdir] & ~endsWith({files.name}, {'.m', '.csv', '.zip', '.png'});
    files = files(valid_idx);
end

num_files = length(files);
fprintf('Found %d tail polar files. Parsing data...\n', num_files);

tail_results = struct();
valid_count = 0;

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
        
        if line_str(1) == '0'
            vals = sscanf(line_str, '%f');
            if length(vals) >= 10 && abs(vals(1)) < 1e-4
                a_val  = vals(2);
                cl_val = vals(5);
                cd_val = vals(6);
                cm_val = vals(10);
                
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
    
    % Extract geometry: St_ARt_lt
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
    
    % Physical Dimensions
    if ~isnan(St_val) && ~isnan(AR_val) && ~isnan(lt_val)
        St_m2  = St_val * 1e-6;
        lt_m   = lt_val * 1e-3;
        Vh_val = (St_m2 * lt_m) / (Sw * c_mac);
        bt_val = sqrt(AR_val * St_val); % Tail span in mm
        ct_val = St_val / bt_val;        % Tail chord in mm
    else
        Vh_val = NaN; bt_val = NaN; ct_val = NaN;
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
    tail_results(valid_count).bt           = bt_val;
    tail_results(valid_count).ct           = ct_val;
    tail_results(valid_count).Vh           = Vh_val;
    tail_results(valid_count).Max_LD       = max_ld;
    tail_results(valid_count).Alpha_Opt    = alpha(idx_ld);
    tail_results(valid_count).CD_min       = cd_min;
    tail_results(valid_count).Max_PF       = max_pf;
    tail_results(valid_count).Cm_alpha     = cm_alpha;
    tail_results(valid_count).alpha        = alpha;
    tail_results(valid_count).CL_CD        = CL_CD;
    tail_results(valid_count).Cm           = Cm;
end

%% Scoring Formulation (Aerodynamics + Stability + Manufacturing Chord)
all_ld   = [tail_results.Max_LD];
all_pf   = [tail_results.Max_PF];
all_cm_a = [tail_results.Cm_alpha];
all_ct   = [tail_results.ct];

norm_ld = (all_ld - min(all_ld)) ./ (max(all_ld) - min(all_ld) + eps);
norm_pf = (all_pf - min(all_pf)) ./ (max(all_pf) - min(all_pf) + eps);
stiff_val = -all_cm_a;
norm_cm = (stiff_val - min(stiff_val)) ./ (max(stiff_val) - min(stiff_val) + eps);
% Wider chord scores higher for buildability and low-Re health
norm_ct = (all_ct - min(all_ct)) ./ (max(all_ct) - min(all_ct) + eps);

for i = 1:valid_count
    % 45% Glide Range + 25% Hang Time + 15% Pitch Stiffness + 15% Manufacturing Chord Depth
    base_score = (0.45 * norm_ld(i) + 0.25 * norm_pf(i) + 0.15 * norm_cm(i) + 0.15 * norm_ct(i)) * 100;
    
    vh_curr = tail_results(i).Vh;
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

%% Display Leaderboard with Dimensions
fprintf('\n================================================================================================================================================\n');
fprintf('%-5s %-28s %-6s %-7s %-7s %-9s %-8s %-16s %-7s %-16s\n', ...
    'Rank', 'Configuration', 'Vh', 'Span(mm)', 'Chord(mm)', 'Max L/D', 'CD_min', 'Cm_alpha(deg^-1)', 'Score', 'Status');
fprintf('================================================================================================================================================\n');
for r = 1:min(20, valid_count)
    fprintf('%-5d %-28s %-6.3f %-9.1f %-9.1f %-9.2f %-8.4f %-16.4f %-7.1f %-16s\n', ...
        r, sorted(r).ConfigName(1:min(28, length(sorted(r).ConfigName))), ...
        sorted(r).Vh, sorted(r).bt, sorted(r).ct, sorted(r).Max_LD, sorted(r).CD_min, ...
        sorted(r).Cm_alpha, sorted(r).Score, sorted(r).Status);
end
fprintf('================================================================================================================================================\n\n');

w = sorted(1);
fprintf('>>> DFM-OPTIMIZED WINNING TAIL VERDICT <<<\n');
fprintf('  Configuration Name    : %s\n', w.ConfigName);
fprintf('  Tail Volume (Vh)      : %.3f\n', w.Vh);
fprintf('  Planform Area (St)    : %.1f mm^2\n', w.St);
fprintf('  Tail Span (bt)        : %.1f mm\n', w.bt);
fprintf('  Tail Chord (ct)       : %.1f mm (Practical Build Thickness)\n', w.ct);
fprintf('  Moment Arm (lt)       : %.1f mm\n', w.lt);
fprintf('  Finite Glide (L/D)    : %.2f\n', w.Max_LD);
fprintf('  Pitch Slope (Cm_a)    : %.4f deg^-1\n', w.Cm_alpha);
fprintf('  Composite DFM Score   : %.1f / 100\n\n', w.Score);

%% Export CSV
csv_name = 'TAIL_ANALYSIS_DFM_RANKED.csv';
fid = fopen(csv_name, 'w');
if fid ~= -1
    fprintf(fid, 'Rank,ConfigName,St_mm2,ARt,lt_mm,bt_mm,ct_mm,Vh,Max_LD,Alpha_Opt,CD_min,Cm_alpha_deg,Score,Status\n');
    for r = 1:valid_count
        it = sorted(r);
        fprintf(fid, '%d,%s,%.1f,%.2f,%.1f,%.1f,%.1f,%.3f,%.2f,%.2f,%.4f,%.4f,%.1f,%s\n', ...
            r, it.ConfigName, it.St, it.ARt, it.lt, it.bt, it.ct, it.Vh, it.Max_LD, it.Alpha_Opt, ...
            it.CD_min, it.Cm_alpha, it.Score, it.Status);
    end
    fclose(fid);
    fprintf('[OK] Full ranked matrix exported to "%s"\n', csv_name);
end