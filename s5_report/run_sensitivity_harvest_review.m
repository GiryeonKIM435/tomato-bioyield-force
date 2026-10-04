function results = run_sensitivity_harvest_review(cfg, opts)
%RUN_SENSITIVITY_HARVEST_REVIEW Refit 4.1–4.4 on each harvest-day cohort
%
% Day subsets use the existing prepare + arm runner. Pooled-model errors split
% by day are assembled in collect_review_sensitivity_tables (no extra replay).

if nargin < 1 || isempty(cfg)
    cfg = PaperStudyConfig();
end
if nargin < 2 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, "reuseExisting")
    opts.reuseExisting = true;
end
if ~isfield(opts, "skipIfExists")
    opts.skipIfExists = true;
end
if ~isfield(opts, "forceRecompute")
    opts.forceRecompute = false;
end

days = ["0423", "0429"];
results = struct();
results.days = days;

for i = 1:numel(days)
    day = days(i);
    fprintf("=== sensitivity harvest: prepare %s ===\n", day);
    run_prepare_jeffreys_harvest_day_subset(cfg, struct( ...
        "day", day, ...
        "skipIfExists", opts.skipIfExists, ...
        "forceRecompute", opts.forceRecompute));

    sensCfg = applyJeffreysHarvestDayPaths(PaperStudyConfig(), day);
    sensCfg.deploy.krVariant = "chord";
    sensCfg = syncPaperKrVariant(sensCfg);
    sensCfg.q7.methodTypes = "force_abs";
    fprintf("=== sensitivity harvest: analyses %s ===\n", sensCfg.analysis.primaryAnalysisTag);
    results.(char("d" + day)) = run_sensitivity_arm_analyses(sensCfg, struct( ...
        "reuseExisting", opts.reuseExisting, ...
        "writeFigures", false));
end
end
