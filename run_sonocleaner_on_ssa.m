%% run_sonocleaner_on_ssa.m
% Point this script at a .ssa file, choose a trace, and run sonocleaner_auto.
%
% REQUIREMENT:
%   sonocleaner_auto.m must be on your MATLAB path or in the same folder.
%
% USAGE:
%   - Set ssaFile below
%   - Either set traceLabelWanted OR traceColumnWanted
%   - Run the script

clear; clc;

%% ---------------- USER SETTINGS ----------------
ssaFile = "tutorial.ssa";     % <-- path to .ssa file

% Choose ONE of these:
traceLabelWanted  = "TRX02:03";   % set to "" if selecting by column number
traceColumnWanted = [];           % e.g. 2 means second data column after time

% Automatic cleaning parameters
maxNoise = 0.02;
levelShiftTolerance = 0.08;
samplingRadius = 4;
verbose = true;
%% -----------------------------------------------

% Parse the SSA file
ssa = parse_ssa_file(ssaFile);

fprintf('Loaded file: %s\n', ssaFile);
fprintf('Sample interval: %.9f s\n', ssa.sampleTimeInterval);
fprintf('Rows: %d\n', size(ssa.data,1));
fprintf('Trace columns found:\n');
for k = 1:numel(ssa.traceLabels)
    fprintf('  %2d : %s\n', k, ssa.traceLabels{k});
end
fprintf('\n');

% Select a trace
[traceVec, traceLabel, traceIndex] = select_trace(ssa, traceLabelWanted, traceColumnWanted);

fprintf('Selected trace %d: %s\n', traceIndex, traceLabel);

% Run automatic cleaner
out = sonocleaner_auto(traceVec, ...
    'MaxNoise', maxNoise, ...
    'LevelShiftTolerance', levelShiftTolerance, ...
    'SamplingRadius', samplingRadius, ...
    'Verbose', verbose);

% Time vector
if ~isempty(ssa.time)
    t = ssa.time;
else
    t = (0:numel(traceVec)-1)' * ssa.sampleTimeInterval;
end

% Plot raw vs cleaned
figure;
hRaw = plot(t, out.x_raw, 'DisplayName', 'Raw'); hold on;
hClean = plot(t, out.x_clean, 'LineWidth', 1.2, 'DisplayName', 'Cleaned');
xlabel('Time (s)');
ylabel('Distance');
title(sprintf('SonoCleaner auto: %s', traceLabel), 'Interpreter', 'none');
legend([hRaw, hClean], {'Raw', 'Cleaned'}, 'Location', 'best');
grid on;

% Plot derivatives and detected shift segments
figure;

subplot(3,1,1);
hRaw1 = plot(t, out.x_raw, 'DisplayName', 'Raw'); hold on;
hClean1 = plot(t, out.x_clean, 'LineWidth', 1.2, 'DisplayName', 'Cleaned');
ylabel('Distance');
title(sprintf('Trace: %s', traceLabel), 'Interpreter', 'none');
legend([hRaw1, hClean1], {'Raw', 'Cleaned'}, 'Location', 'best');
grid on;

subplot(3,1,2);
td = t(1:end-1);
hDiffRaw = plot(td, out.d_raw, 'DisplayName', 'Raw diff'); hold on;
hDiffClean = plot(td, out.d_clean, 'LineWidth', 1.2, 'DisplayName', 'Clean diff');

if ~isempty(out.labelled_segments)
    xl = xline(td(out.labelled_segments), '--');
    for k = 1:numel(xl)
        xl(k).HandleVisibility = 'off';   % hide vertical lines from legend
    end
end

ylabel('\Deltax');
legend([hDiffRaw, hDiffClean], {'Raw diff', 'Clean diff'}, 'Location', 'best');
grid on;

subplot(3,1,3);
tdd = t(1:end-2);
plot(tdd, out.dd_raw, 'HandleVisibility', 'off');
ylabel('\Delta^2x');
xlabel('Time (s)');
grid on;

fprintf('\nDone.\n');
fprintf('Labelled segments: %d\n', numel(out.labelled_segments));
fprintf('Matched groups   : %d\n', numel(out.matches));

%% ================= LOCAL FUNCTIONS =================

function ssa = parse_ssa_file(filename)
    % Parses a Sonometrics .ssa file of the style used by SonoCleaner.
    %
    % Returns struct with fields:
    %   .sampleTimeInterval
    %   .headersRaw
    %   .traceLabels
    %   .units
    %   .time
    %   .data        (numeric matrix of trace columns only, excluding time)
    %   .allNumeric   (full numeric table including time column)

    txt = fileread(filename);
    lines = regexp(txt, '\r\n|\n|\r', 'split');
    lines = lines(:);

    sampleTimeInterval = [];
    beginIdx = [];
    endIdx = [];

    for i = 1:numel(lines)
        line = strtrim(lines{i});

        if startsWith(line, 'SAMPLE TIME INTERVAL:', 'IgnoreCase', true)
            val = strtrim(extractAfter(line, 'SAMPLE TIME INTERVAL:'));
            sampleTimeInterval = str2double(val);
        elseif strcmpi(line, 'BEGIN DATA:') || strcmpi(line, 'BEGIN DATA')
            beginIdx = i;
        elseif strcmpi(line, 'END DATA') || strcmpi(line, 'END DATA:')
            endIdx = i;
            break;
        end
    end

    if isempty(beginIdx)
        error('Could not find BEGIN DATA in %s', filename);
    end
    if isempty(endIdx)
        error('Could not find END DATA in %s', filename);
    end
    if isempty(sampleTimeInterval) || isnan(sampleTimeInterval)
        warning('Could not parse SAMPLE TIME INTERVAL. Falling back to empty.');
        sampleTimeInterval = [];
    end

    % In SSA, typically:
    % line beginIdx+1 = headers
    % line beginIdx+2 = units
    % line beginIdx+3:endIdx-1 = numeric data
    if beginIdx + 2 > endIdx - 1
        error('Data block appears malformed.');
    end

    headerLine = lines{beginIdx+1};
    unitLine   = lines{beginIdx+2};

    headersRaw = split_tab_preserve_nonempty(headerLine);
    units      = split_tab_preserve_nonempty(unitLine);

    dataLines = lines(beginIdx+3:endIdx-1);
    dataLines = dataLines(~cellfun(@(c) isempty(strtrim(c)), dataLines));

    numericRows = cell(numel(dataLines), 1);
    maxCols = 0;

    for i = 1:numel(dataLines)
        toks = split_tab_preserve_nonempty(dataLines{i});
        vals = nan(1, numel(toks));
        for j = 1:numel(toks)
            vals(j) = str2double(strtrim(toks{j}));
        end
        numericRows{i} = vals;
        maxCols = max(maxCols, numel(vals));
    end

    allNumeric = nan(numel(numericRows), maxCols);
    for i = 1:numel(numericRows)
        row = numericRows{i};
        allNumeric(i,1:numel(row)) = row;
    end

    if size(allNumeric,2) < 2
        error('No trace columns found in SSA file.');
    end

    time = allNumeric(:,1);
    data = allNumeric(:,2:end);

    % Strip any trailing junk labels like "Delineators"
    traceLabels = headersRaw(2:min(numel(headersRaw), size(allNumeric,2)));
    traceLabels = cellfun(@strtrim, traceLabels, 'UniformOutput', false);

    % Keep only as many labels as there are numeric trace columns
    if numel(traceLabels) > size(data,2)
        traceLabels = traceLabels(1:size(data,2));
    elseif numel(traceLabels) < size(data,2)
        % pad with generic labels if needed
        for k = numel(traceLabels)+1:size(data,2)
            traceLabels{k} = sprintf('Trace_%d', k);
        end
    end

    ssa = struct();
    ssa.sampleTimeInterval = sampleTimeInterval;
    ssa.headersRaw = headersRaw;
    ssa.traceLabels = traceLabels;
    ssa.units = units;
    ssa.time = time;
    ssa.data = data;
    ssa.allNumeric = allNumeric;
end

function [traceVec, traceLabel, traceIndex] = select_trace(ssa, traceLabelWanted, traceColumnWanted)
    traceLabelWanted = string(traceLabelWanted);

    if strlength(traceLabelWanted) > 0
        idx = find(strcmpi(strtrim(ssa.traceLabels), strtrim(traceLabelWanted)), 1);
        if isempty(idx)
            error('Trace label "%s" not found.', traceLabelWanted);
        end
        traceIndex = idx;
    elseif ~isempty(traceColumnWanted)
        traceIndex = traceColumnWanted;
        if traceIndex < 1 || traceIndex > size(ssa.data,2)
            error('Trace column %d is out of range.', traceIndex);
        end
    else
        error('Set either traceLabelWanted or traceColumnWanted.');
    end

    traceVec = ssa.data(:, traceIndex);
    traceLabel = ssa.traceLabels{traceIndex};
end

function toks = split_tab_preserve_nonempty(line)
    % SSA files are tab-delimited, often with repeated tabs.
    toks = regexp(line, '\t', 'split');
    toks = toks(~cellfun(@isempty, toks));
end