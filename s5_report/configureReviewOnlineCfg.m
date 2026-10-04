function cfg = configureReviewOnlineCfg(tag, opts)
%configureReviewOnlineCfg Sequential-replay cfg for one review sensitivity arm
%
% Output tag is separate from the manuscript tag. Cohort mats stay on the
% shipped Jeffreys IQR paths. Trajectory reuse is opt-in via opts.reuseTrajTag.

if nargin < 2 || isempty(opts)
    opts = struct();
end

cfg = PaperStudyConfig();
cfg.deploy.krVariant = "chord";
cfg = syncPaperKrVariant(cfg);
cfg.q7.methodTypes = "force_abs";

tag = char(string(tag));
cfg.q7.analysisTag = tag;
cfg.analysis.primaryAnalysisTag = tag;
cfg.paper.q3AnalysisTag = tag;

if isfield(opts, "minBandPoints") && isfinite(opts.minBandPoints)
    cfg.deploy.minBandPointsForKr = opts.minBandPoints;
end
if isfield(opts, "quantileP") && isfinite(opts.quantileP)
    cfg.q7.quantileP = opts.quantileP;
end
if isfield(opts, "reuseTrajTag") && strlength(string(opts.reuseTrajTag)) > 0
    cfg.q7.reuseQ3TrajTag = char(string(opts.reuseTrajTag));
    cfg.cache.cohortAnalysisTag = char(string(opts.reuseTrajTag));
else
    cfg.cache.cohortAnalysisTag = tag;
end
end
