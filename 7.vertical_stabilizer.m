% =========================================================================
% FLOW5 VERTICAL TAIL DFM & AERODYNAMIC SCORING SCRIPT (WITH DRAG)
% Department of Mechanical Engineering, University of Moratuwa
% =========================================================================
clear; clc; close all;

%% 1. Path Setup & Folder Ingestion
% Target directory where your flow5 polar txt files are located
folderPath = fullfile('Fin_analysis_polars', 'Fin_analysis_polars');

% Fallback verification: if nested folder does not exist, check primary folder
if ~exist(folderPath, 'dir')
    if exist('Fin_analysis_polars', 'dir')
        folderPath = 'Fin_analysis_polars';
    else
        error('Cannot find folder "%s". Please verify the path.', folderPath);
    end
end

filePattern = fullfile(folderPath, '*.txt'); 
files = dir(filePattern);

numFiles = length(files);
if numFiles == 0
    error('No .txt files found in directory: %s', folderPath);
end

Model_Name    = strings(numFiles, 1);
File_Name     = strings(numFiles, 1); 
Cn_beta_Slope = NaN(numFiles, 1); 
Tail_Area     = NaN(numFiles, 1); 
Mean_CD       = NaN(numFiles, 1);

fprintf('Reading %d polar data files from "%s"...\n', numFiles, folderPath);

%% 2. Data Extraction Loop
for i = 1:numFiles
    fileName     = files(i).name;
    fullFilePath = fullfile(folderPath, fileName);
    
    File_Name(i) = string(fileName);
    [~, name, ~] = fileparts(fileName);
    Model_Name(i) = string(name);
    
    % Parse area from file prefix (e.g., "0_04_..." -> 0.04)
    nameParts = split(Model_Name(i), '_');
    if length(nameParts) >= 2
        Tail_Area(i) = str2double(nameParts(1) + "." + nameParts(2));
    end
    
    fid = fopen(fullFilePath, 'rt');
    if fid == -1, continue; end
    
    betaCol = -1; cnCol = -1; cdCol = -1; headerFound = false;
    
    % Search for flow5 header definition
    while ~feof(fid)
        tline = strtrim(fgetl(fid));
        if isempty(tline), continue; end
        
        if contains(tline, 'Ctrl') && contains(tline, 'Cn')
            headerFound = true;
            cleanLine = regexprep(tline, '\([^)]*\)', ''); 
            tokens = strsplit(strtrim(cleanLine));
            
            for k = 1:length(tokens)
                tok = lower(tokens{k});
                if tok == "b" || tok == "beta" || contains(tok, char(946))
                    betaCol = k;
                elseif strcmp(tokens{k}, 'Cn')
                    cnCol = k;
                elseif strcmp(tokens{k}, 'CD') || strcmp(tokens{k}, 'Cd')
                    cdCol = k;
                end
            end
            
            % Default flow5 column indices if header string split shifts
            if betaCol == -1, betaCol = 3; end
            if cdCol == -1,   cdCol = 6;   end
            if cnCol == -1,   cnCol = 14;  end
            break;
        end
    end
    
    % Extract numeric operating points
    beta_vals = []; cn_vals = []; cd_vals = [];
    if headerFound
        while ~feof(fid)
            tline = strtrim(fgetl(fid));
            if isempty(tline), continue; end
            
            dataTokens = strsplit(tline);
            maxColRequired = max([betaCol, cnCol, cdCol]);
            
            if length(dataTokens) >= maxColRequired
                b_num = str2double(dataTokens{betaCol});
                c_num = str2double(dataTokens{cnCol});
                d_num = str2double(dataTokens{cdCol});
                
                if ~isnan(b_num) && ~isnan(c_num) && ~isnan(d_num)
                    beta_vals(end+1) = b_num; %#ok<AGROW>
                    cn_vals(end+1)   = c_num; %#ok<AGROW>
                    cd_vals(end+1)   = d_num; %#ok<AGROW>
                end
            end
        end
    end
    fclose(fid);
    
    % Calculate directional derivative (dCn/dBeta) and average drag
    if length(beta_vals) > 1 && length(cn_vals) > 1
        p = polyfit(beta_vals, cn_vals, 1);
        Cn_beta_Slope(i) = p(1);       % Flow5 raw slope (deg^-1)
        Mean_CD(i)       = mean(cd_vals); 
    end
end

resultsTable = table(Model_Name, File_Name, Tail_Area, Mean_CD, Cn_beta_Slope);
resultsTable(isnan(resultsTable.Cn_beta_Slope), :) = []; 

%% 3. Multi-Attribute Scoring & Ranking
if height(resultsTable) > 0
    % 1. Stability Score (50% Weight):
    % In flow5 body axis, more negative slope = stronger weathercock restoring torque.
    % Most negative slope earns 100 points; least negative earns 0 points.
    min_Cn = min(resultsTable.Cn_beta_Slope);
    max_Cn = max(resultsTable.Cn_beta_Slope);
    if abs(max_Cn - min_Cn) > eps
        resultsTable.Stability_Score = ((max_Cn - resultsTable.Cn_beta_Slope) ./ (max_Cn - min_Cn)) * 100;
    else
        resultsTable.Stability_Score = repmat(100, height(resultsTable), 1);
    end
    
    % 2. DFM / Area Score (30% Weight):
    % Smaller fin area reduces structural weight on tail boom.
    min_Area = min(resultsTable.Tail_Area);
    max_Area = max(resultsTable.Tail_Area);
    if abs(max_Area - min_Area) > eps
        resultsTable.DFM_Score = ((max_Area - resultsTable.Tail_Area) ./ (max_Area - min_Area)) * 100;
    else
        resultsTable.DFM_Score = repmat(100, height(resultsTable), 1);
    end

    % 3. Drag Score (20% Weight):
    % Lower average drag floor earns higher score.
    min_CD = min(resultsTable.Mean_CD);
    max_CD = max(resultsTable.Mean_CD);
    if abs(max_CD - min_CD) > eps
        resultsTable.Drag_Score = ((max_CD - resultsTable.Mean_CD) ./ (max_CD - min_CD)) * 100;
    else
        resultsTable.Drag_Score = repmat(100, height(resultsTable), 1);
    end
    
    % Composite Score Formulation (50% Stability + 30% DFM + 20% Drag)
    resultsTable.Total_Marks = (0.50 * resultsTable.Stability_Score) + ...
                               (0.30 * resultsTable.DFM_Score) + ...
                               (0.20 * resultsTable.Drag_Score);
    
    % Sort by Total Marks descending
    rankedTable = sortrows(resultsTable, 'Total_Marks', 'descend');

    fprintf('\n================== TOP 10 VERTICAL TAIL CANDIDATES ==================\n');
    disp(head(rankedTable(:, {'Model_Name', 'Stability_Score', 'DFM_Score', 'Drag_Score', 'Total_Marks'}), 10));

    % Export CSV Table
    outputFilename = 'Ranked_Vertical_Tails_Complete.csv';
    writetable(rankedTable(:, {'Model_Name', 'Tail_Area', 'Mean_CD', 'Cn_beta_Slope', ...
                               'Stability_Score', 'DFM_Score', 'Drag_Score', 'Total_Marks'}), outputFilename);
    fprintf('\nFull ranked results exported successfully to: "%s"\n', outputFilename);

    % Display Winner
    bestModel = rankedTable(1, :);
    fprintf('\n********************************************************\n');
    fprintf('   WINNING CONFIGURATION : %s\n', bestModel.Model_Name);
    fprintf('   Total Composite Score : %.2f / 100\n', bestModel.Total_Marks);
    fprintf('   Breakdown             : [Stability: %.1f | DFM: %.1f | Drag: %.1f]\n', ...
        bestModel.Stability_Score, bestModel.DFM_Score, bestModel.Drag_Score);
    fprintf('   dCn/dBeta Slope       : %.6f deg^-1\n', bestModel.Cn_beta_Slope);
    fprintf('   Mean CD               : %.6f\n', bestModel.Mean_CD);
    fprintf('********************************************************\n');
else
    fprintf('\nError: No valid data could be processed from polar files.\n');
end