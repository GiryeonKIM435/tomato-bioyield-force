function results = run_sensitivity_kupdate(cfg, opts)
%RUN_SENSITIVITY_KUPDATE Sequential replay at alternate k-update point gates
%
% Manuscript gate is 5 points (jeffreys_bi_iqr15) and is not recomputed.
% Default arms: 2 and 20 points -> sens_kupdate_n2 / sens_kupdate_n20.

if nargin < 1 || isempty(cfg)
    cfg = PaperStudyConfig();
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, "reuseExisting")
    opts.reuseExisting = true;
end
if ~isfield(opts, "pointCounts") || isempty(opts.pointCounts)
    opts.pointCounts = [2, 20];
end

pointCounts = double(opts.pointCounts(:));
results = cell(numel(pointCounts), 1);
for i = 1:numel(pointCounts)
    nPts = pointCounts(i);
    tag = sprintf("sens_kupdate_n%d", nPts);
    cfgRun = configureReviewOnlineCfg(tag, struct("minBandPoints", nPts));
    fprintf("=== sensitivity k-update: minBandPoints=%d -> %s ===\n", nPts, tag);
    results{i} = run_online_evaluation(cfgRun, struct( ...
        "useOutlierFilter", false, ...
        "analysisTag", tag, ...
        "reuseExisting", opts.reuseExisting, ...
        "writeFigures", false));
    results{i}.minBandPoints = nPts;
    results{i}.manuscriptTag = "jeffreys_bi_iqr15";
end
end
