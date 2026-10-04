function results = run_sensitivity_quantile(cfg, opts)
%RUN_SENSITIVITY_QUANTILE Sequential replay at relative-overestimation p=0.90 and 0.99
%
% p=0.95 is the manuscript run (jeffreys_bi_iqr15) and is not recomputed.
% Stiffness trajectories are reused from that tag. Stop thresholds are not.

if nargin < 1 || isempty(cfg)
    cfg = PaperStudyConfig();
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, "reuseExisting")
    opts.reuseExisting = true;
end

manuscriptTag = "jeffreys_bi_iqr15";
levels = [0.90, 0.99];
tags = ["sens_quantile_p90", "sens_quantile_p99"];
results = cell(numel(levels), 1);

for i = 1:numel(levels)
    tag = char(tags(i));
    cfgRun = configureReviewOnlineCfg(tag, struct( ...
        "quantileP", levels(i), ...
        "minBandPoints", 5, ...
        "reuseTrajTag", manuscriptTag));
    fprintf("=== sensitivity quantile: p=%.2f -> %s (traj=%s) ===\n", ...
        levels(i), tag, manuscriptTag);
    results{i} = run_online_evaluation(cfgRun, struct( ...
        "useOutlierFilter", false, ...
        "analysisTag", tag, ...
        "reuseExisting", opts.reuseExisting, ...
        "writeFigures", false));
    results{i}.quantileP = levels(i);
end
end
