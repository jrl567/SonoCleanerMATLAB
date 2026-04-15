function out = sonocleaner_auto(x, varargin)
% SONOCLEANER_AUTO
% MATLAB implementation made as faithful as reasonably possible to the
% automatic procedure in SonoCleaner's Haskell source.
%
% INPUT
%   x : Nx1 or 1xN trace
%
% NAME-VALUE PARAMETERS
%   'MaxNoise'             : noise threshold (default 0.02)
%   'LevelShiftTolerance'  : level-shift tolerance for labelling (default 0.08)
%   'SamplingRadius'       : local slope estimation radius (default 4)
%   'MaxNonLevelShifts'    : max low-slope segments inside labelled region (default 4)
%   'GroupCountLimit'      : max accepted zero-sum group size (default 8)
%   'SearchGroupCountLimit': max search size (default 16)
%   'InterpolationLimit'   : short-span interpolation limit (default 3)
%   'Verbose'              : true/false (default false)
%
% OUTPUT
%   out : struct with fields
%       .x_raw
%       .x_clean
%       .d_raw
%       .d_clean
%       .dd_raw
%       .labelled_segments
%       .slope_estimates
%       .slope_errors
%       .matches
%       .updates
%
% NOTES
% - This only implements the automated procedure, not the manual tools.
% - Remaining small differences from Haskell are mostly due to MATLAB vs
%   Haskell library behavior (especially signrank/Wilcoxon details).

    p = inputParser;
    addParameter(p, 'MaxNoise', 0.02, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
    addParameter(p, 'LevelShiftTolerance', 0.08, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
    addParameter(p, 'SamplingRadius', 4, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
    addParameter(p, 'MaxNonLevelShifts', 4, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
    addParameter(p, 'GroupCountLimit', 8, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
    addParameter(p, 'SearchGroupCountLimit', 16, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
    addParameter(p, 'InterpolationLimit', 3, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
    addParameter(p, 'Verbose', false, @(v)islogical(v)&&isscalar(v));
    parse(p, varargin{:});
    cfg = p.Results;

    x = x(:);
    n = numel(x);
    if n < 3
        error('Trace must have at least 3 samples.');
    end

    % TraceState-style core series
    d = diff(x);      % diffSeries
    dd = diff(d);     % diff2Series inner vector

    % 1) Label level-shifts
    labelled = label_level_shifts(d, dd, cfg.MaxNoise, cfg.LevelShiftTolerance, cfg.MaxNonLevelShifts);

    % 2) Estimate slopes and slope-errors exactly in the same spirit as Matching.hs
    slope_est = zeros(size(d));
    slope_err = nan(size(d));
    badSegments = labelled(:).';  % no manual modifiedSegments here, so badSegments = labelled only

    for i = labelled(:).'
        slope_est(i) = estimate_local_slope(d, badSegments, cfg.SamplingRadius, i);
        slope_err(i) = d(i) - slope_est(i);
    end

    % 3) Match zero-sum groups with Haskell-style candidate span search
    matches = match_level_shifts_haskell_style( ...
        x, labelled, slope_est, slope_err, ...
        cfg.MaxNoise, cfg.GroupCountLimit, cfg.SearchGroupCountLimit, cfg.InterpolationLimit);

    % 4) Apply all matches, equivalent to applyMatches matches progression
    %    with progression = all groups.
    d_clean = d;
    updates = zeros(0,2);

    for k = 1:numel(matches)
        update_pairs = matches(k).updates;
        if isempty(update_pairs)
            continue;
        end
        d_clean(update_pairs(:,1)) = update_pairs(:,2);
        updates = [updates; update_pairs];
    end

    % 5) Reconstruct the cleaned series from the cleaned diff series
    x_clean = reconstruct_from_diff(x(1), d_clean);

    out = struct();
    out.x_raw = x;
    out.x_clean = x_clean;
    out.d_raw = d;
    out.d_clean = d_clean;
    out.dd_raw = dd;
    out.labelled_segments = labelled(:);
    out.slope_estimates = slope_est;
    out.slope_errors = slope_err;
    out.matches = matches;
    out.updates = updates;

    if cfg.Verbose
        fprintf('Labelled segments: %d\n', numel(labelled));
        fprintf('Matched groups   : %d\n', numel(matches));
    end
end

% =========================================================================
% LABELING: close port of Labelling.hs
% =========================================================================
function labelled = label_level_shifts(d, dd, maxNoise, levelShiftTolerance, maxNonLevelShifts)
% Mirrors labelLevelShifts' in Labelling.hs.

    slopeLimit = levelShiftTolerance - maxNoise;
    curveLimit = levelShiftTolerance - maxNoise;
    matchLimit = 2 * maxNoise;

    if isempty(dd)
        labelled = [];
        return;
    end

    % Haskell uses ivExtend1 1 and ivExtend2 1 before searching.
    % For dv this is equivalent to padding the first-difference series by
    % repeating the end values once at each edge.
    d_ext = [d(1); d(:); d(end)];
    dd_ext = diff(d_ext);  % extended second-difference series

    high_curv = find(abs(dd_ext) > curveLimit);  % analogous to ivFindIndices2
    matched_intervals = zeros(0,2);

    k = 1;
    while k <= numel(high_curv)
        i0 = high_curv(k);
        is = high_curv((k+1):end);
        found = false;

        for t = 1:numel(is)
            i1 = is(t);

            % Candidate IndexInterval on the extended dd index set.
            % innerInterval = iiShrink . iiUndiff
            % For interval (i0, i1), interior slope indices are (i0+1):i1.
            s0 = i0 + 1;
            s1 = i1;

            if s0 > s1 || s1 > numel(d_ext)
                continue;
            end

            % countNonLevelShifts = ivCount ((< slopeLimit) . abs) . (`ivSlice` dv)
            nonLevelCount = sum(abs(d_ext(s0:s1)) < slopeLimit);

            % V.takeWhile means that once the count exceeds the limit, later
            % wider intervals are not considered for this i0.
            if nonLevelCount > maxNonLevelShifts
                break;
            end

            % V.map (uncurry subtract . iiIndex dv . iiUndiff)
            % iiUndiff maps dd interval (i0,i1) to dv endpoint interval (i0, i1+1),
            % so total change in slope is dv(i1+1) - dv(i0).
            if (i1 + 1) > numel(d_ext)
                continue;
            end
            totalSlopeChange = d_ext(i1 + 1) - d_ext(i0);

            % First interval passing this test is chosen.
            if abs(totalSlopeChange) <= matchLimit
                matched_intervals(end+1, :) = [i0, i1];
                k = k + t + 1;  % V.drop (n+1) is
                found = true;
                break;
            end
        end

        if ~found
            k = k + 1;
        end
    end

    % Concatenate inner intervals and undo the 1-sample extension shift.
    labelled_mask_ext = false(numel(d_ext), 1);
    for r = 1:size(matched_intervals, 1)
        i0 = matched_intervals(r,1);
        i1 = matched_intervals(r,2);
        s0 = i0 + 1;
        s1 = i1;
        labelled_mask_ext(s0:s1) = true;
    end

    % Original d corresponds to d_ext(2:end-1)
    labelled = find(labelled_mask_ext(2:end-1));
end

% =========================================================================
% SLOPE ESTIMATION: close port of Slope.hs
% =========================================================================
function s = estimate_local_slope(d, badSegments, radius, i)
% Mirrors estimateSlope in Slope.hs:
%   - surrounding slopes within radius
%   - exclude bad segments
%   - nonzero subset used for Wilcoxon
%   - if <5 nonzero, return 0
%   - if p<0.10, return median of all surrounding slopes
%   - otherwise return 0

    lo = max(1, i - radius);
    hi = min(numel(d), i + radius);
    idx = lo:hi;
    idx = idx(~ismember(idx, badSegments));

    slopes = d(idx);
    nonzero = slopes(slopes ~= 0);

    if numel(nonzero) < 5
        s = 0;
        return;
    end

    try
        p = signrank(nonzero, 0, 'tail', 'both');
    catch
        % Fallback if signrank chokes on ties.
        p = approximate_signrank_pvalue(nonzero);
    end

    if p < 0.10
        s = median(slopes);
    else
        s = 0;
    end
end

function p = approximate_signrank_pvalue(x)
% Crude fallback only. Main path should be MATLAB signrank.
    x = x(:);
    x = x(x ~= 0);
    n = numel(x);
    if n < 5
        p = 1;
        return;
    end
    pos = sum(x > 0);
    % Very rough two-sided sign test fallback
    p_one = binocdf(min(pos, n-pos), n, 0.5);
    p = min(1, 2 * p_one);
end

% =========================================================================
% MATCHING: sortof close port of Matching.hs
% =========================================================================
function matches = match_level_shifts_haskell_style(x, labelled, slope_est, slope_err, noiseTh, groupCountLimit, searchGroupCountLimit, interpolationLimit)
% Follows Matching.hs mostly:
%   - initialize heap with consecutive distances
%   - pop smallest (span, startID)
%   - search for exact target span from that start
%   - if overshoot, reschedule with next larger span
%   - if zero-sum found, remove those level-shifts from the chain
%   - continue until heap empty

    if isempty(labelled)
        matches = struct('positions', {}, 'span', {}, 'errorSum', {}, ...
                         'method', {}, 'updates', {}, 'visitedIDs', {});
        return;
    end

    % "IndexedChain" equivalent
    chain.pos = labelled(:).';
    chain.err = slope_err(labelled(:)).';
    chain.exists = true(size(chain.pos));
    chain.next = [2:numel(chain.pos), 0];   % 0 means nullIndex
    if isempty(chain.next)
        chain.next = 0;
    end

    % Heap seed candidates: distances between consecutive labelled positions,
    % paired with the ID of the first element in the pair.
    % Haskell uses [0..] in 0-based indexing; here use 1-based IDs.
    if numel(labelled) >= 2
        distances = diff(labelled(:)).';
        heap = [distances(:), (1:numel(distances)).'];
    else
        heap = zeros(0,2);
    end

    matchesFlat = struct('positions', {}, 'span', {}, 'errorSum', {}, ...
                         'method', {}, 'updates', {}, 'visitedIDs', {});

    while ~isempty(heap)
        heap = sortrows(heap, [1 2]);   % min-heap behavior: span first, then earlier startID
        targetSpan = heap(1,1);
        targetID   = heap(1,2);
        heap(1,:) = [];

        [searchResult, visitedIDs, groupPositions, groupErr, nextSpan] = ...
            search_zero_sum(targetSpan, targetID, chain, noiseTh, ...
                            groupCountLimit, searchGroupCountLimit);

        switch searchResult
            case 'NoSolution'
                % do nothing

            case 'NextSpan'
                heap(end+1, :) = [nextSpan, targetID];

            case 'ZeroSum'
                % Remove all matched level-shifts from the chain
                for id = visitedIDs(:).'
                    chain = chain_remove(chain, id);
                end

                rec.positions = groupPositions(:).';
                rec.span = targetSpan;
                rec.errorSum = groupErr;
                rec.visitedIDs = visitedIDs(:).';

                if rec.span <= interpolationLimit
                    rec.method = 'interpolate';
                    rec.updates = interpolate_updates_from_signal(x, rec.positions(1), rec.positions(end));
                else
                    rec.method = 'redistribute';
                    new_slopes = slope_est(rec.positions);
                    new_slopes(1) = new_slopes(1) + rec.errorSum;
                    rec.updates = [rec.positions(:), new_slopes(:)];
                end

                matchesFlat(end+1) = rec;
        end
    end

    % Haskell groups matches by span and those span groups increase monotonically.
    % For convenience, keep output flat but sorted the same way.
    if isempty(matchesFlat)
        matches = matchesFlat;
        return;
    end

    spans = [matchesFlat.span].';
    starts = arrayfun(@(m)m.positions(1), matchesFlat).';
    [~, ord] = sortrows([spans, starts], [1 2]);
    matches = matchesFlat(ord);
end

function [resultType, visitedIDs, groupPositions, groupErr, nextSpan] = ...
    search_zero_sum(spanLimit, startIdx, chain, errLim, groupCountLimit, searchGroupCountLimit)
% Mirrors searchZeroSum + searchZeroSumHelper + foldChainFrom.

    resultType = 'NoSolution';
    visitedIDs = [];
    groupPositions = [];
    groupErr = NaN;
    nextSpan = NaN;

    q = chain_query(chain, startIdx);
    if isempty(q)
        return;
    end

    startPos = q.pos;
    targetPos = startPos + spanLimit;

    accCount = 0;
    accErr = 0;
    accPos = [];

    idx = startIdx;

    while idx ~= 0
        q = chain_query(chain, idx);
        if isempty(q)
            return;
        end

        pos = q.pos;
        err = q.err;

        accCount = accCount + 1;
        accErr = accErr + err;
        accPos(end+1) = pos;
        visitedIDs(end+1) = idx;

        if pos < targetPos
            idx = chain_next(chain, idx);
            continue;
        elseif pos == targetPos
            if accCount > searchGroupCountLimit
                resultType = 'NoSolution';
                return;
            elseif accCount > groupCountLimit
                idx = chain_next(chain, idx);
                continue;
            elseif abs(accErr) <= errLim
                resultType = 'ZeroSum';
                groupPositions = accPos;
                groupErr = accErr;
                return;
            else
                idx = chain_next(chain, idx);
                continue;
            end
        else
            resultType = 'NextSpan';
            nextSpan = pos - startPos;
            return;
        end
    end

    resultType = 'NoSolution';
end

function q = chain_query(chain, idx)
    if idx < 1 || idx > numel(chain.pos) || ~chain.exists(idx)
        q = [];
    else
        q.pos = chain.pos(idx);
        q.err = chain.err(idx);
    end
end

function nextIdx = chain_next(chain, idx)
% "next" in IndexedChain skips removed entries.
    if idx < 1 || idx > numel(chain.next)
        nextIdx = 0;
        return;
    end

    nextIdx = chain.next(idx);
    while nextIdx ~= 0
        if chain.exists(nextIdx)
            return;
        end
        nextIdx = chain.next(nextIdx);
    end
end

function chain = chain_remove(chain, idx)
% Equivalent enough to IC.remove for this usage.
    if idx >= 1 && idx <= numel(chain.exists)
        chain.exists(idx) = false;
    end
end

% =========================================================================
% INTERPOLATION / RECONSTRUCTION
% =========================================================================
function update_pairs = interpolate_updates_from_signal(x, firstSeg, lastSeg)
% Equivalent in intent to:
%   interpolationUpdates v (iiUndiff $ IndexInterval (i, i+groupSpan))
%
% If slope segments firstSeg:lastSeg are corrected, the signal points involved
% are x(firstSeg : lastSeg+1). Endpoints are preserved and the interior is
% linearly interpolated.

    leftPt = firstSeg;
    rightPt = lastSeg + 1;

    if leftPt < 1 || rightPt > numel(x) || rightPt <= leftPt
        update_pairs = zeros(0,2);
        return;
    end

    y_left = x(leftPt);
    y_right = x(rightPt);

    xi = leftPt:rightPt;
    yi = linspace(y_left, y_right, numel(xi)).';

    d_interp = diff(yi);
    seg_idx = firstSeg:lastSeg;

    update_pairs = [seg_idx(:), d_interp(:)];
end

function x = reconstruct_from_diff(x0, d)
% Equivalent to ivUndiff on (offset, diffSeries) with offset = x0.
    x = zeros(numel(d)+1, 1);
    x(1) = x0;
    x(2:end) = x0 + cumsum(d(:));
end