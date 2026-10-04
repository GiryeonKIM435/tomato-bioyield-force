function results = export_review_comparison_tables(cfg)
%EXPORT_REVIEW_COMPARISON_TABLES Three comparison tables for the review checks
%
%   compare_quantile.csv          stopping-margin quantile 95% vs 90% vs 99%
%   compare_kupdate_points.csv    chord k update after 5 vs 2 vs 20 samples
%   compare_harvest_day.csv       all harvest days vs day 1 vs day 2

if nargin < 1 || isempty(cfg)
    cfg = PaperStudyConfig();
end

outDir = fullfile(cfg.out.sensitivity, "review");
if ~isfolder(outDir)
    mkdir(outDir);
end
clearPreviousReviewTables(outDir);

gammaTag = "gamma_" + strrep(sprintf("%.1f", cfg.q7.gammaValues(1)), ".", "p");
lockedKey = "force_s05_w30";

quantileSpecs = { ...
    "95% (manuscript)", 0.95, "jeffreys_bi_iqr15"; ...
    "90%", 0.90, "sens_quantile_p90"; ...
    "99%", 0.99, "sens_quantile_p99"};
qRows = cell(size(quantileSpecs, 1), 1);
for i = 1:size(quantileSpecs, 1)
    qRows{i} = sequentialComparisonRow(cfg, gammaTag, lockedKey, ...
        string(quantileSpecs{i, 1}), quantileSpecs{i, 2}, string(quantileSpecs{i, 3}), ...
        "quantile");
end
quantileTable = vertcat(qRows{:});
quantileTable = movevars(quantileTable, "settingValue", "After", "setting");
quantileTable = renamevars(quantileTable, "settingValue", "quantile");

pointSpecs = { ...
    "5 samples (manuscript)", 5, "jeffreys_bi_iqr15"; ...
    "2 samples", 2, "sens_kupdate_n2"; ...
    "20 samples", 20, "sens_kupdate_n20"};
pRows = cell(size(pointSpecs, 1), 1);
for i = 1:size(pointSpecs, 1)
    pRows{i} = sequentialComparisonRow(cfg, gammaTag, lockedKey, ...
        string(pointSpecs{i, 1}), pointSpecs{i, 2}, string(pointSpecs{i, 3}), ...
        "minBandPoints");
end
pointTable = vertcat(pRows{:});
pointTable = movevars(pointTable, "settingValue", "After", "setting");
pointTable = renamevars(pointTable, "settingValue", "minBandPoints");

harvestTable = harvestComparisonTable(cfg, gammaTag, lockedKey);

qPath = fullfile(outDir, "compare_quantile.csv");
kPath = fullfile(outDir, "compare_kupdate_points.csv");
hPath = fullfile(outDir, "compare_harvest_day.csv");
writetable(quantileTable, qPath);
writetable(pointTable, kPath);
writetable(harvestTable, hPath);

results = struct();
results.quantile = quantileTable;
results.kupdate = pointTable;
results.harvest = harvestTable;
results.quantilePath = qPath;
results.kupdatePath = kPath;
results.harvestPath = hPath;
fprintf("Comparison tables:\n  %s\n  %s\n  %s\n", qPath, kPath, hPath);
end

function row = sequentialComparisonRow(cfg, gammaTag, lockedKey, setting, settingValue, tag, ~)
base = fullfile(cfg.out.q7, tag, gammaTag);
summary = readIfExists(fullfile(base, "q7_design_deploy_summary.csv"));
early = readIfExists(fullfile(base, "table_online_early_stop_by_method.csv"));
alphaTbl = readIfExists(fullfile(base, "q7_design_alpha_by_method.csv"));
if isempty(summary)
    error("export_review_comparison_tables:MissingSeq", ...
        "Sequential summary missing for %s (%s).", setting, tag);
end
summary = summary(string(summary.methodType) == "force_abs", :);
locked = metricsAtKey(summary, early, alphaTbl, lockedKey);
[~, ix] = min(summary.finalUpdateMae);
bestKey = string(summary.krMethodKey(ix));
best = metricsAtKey(summary, early, alphaTbl, bestKey);

row = table(setting, settingValue, tag, lockedKey, ...
    locked.mae, locked.maeSem, locked.r2, locked.bioyield, locked.premature, ...
    locked.safeStopRate, locked.relErr, locked.relErrSem, locked.stopMae, ...
    locked.ep, locked.alpha, locked.alphaMin, locked.alphaMax, ...
    bestKey, string(summary.label(ix)), ...
    best.mae, best.maeSem, best.r2, best.bioyield, best.premature, ...
    best.safeStopRate, best.alpha, ...
    'VariableNames', {'setting', 'settingValue', 'analysisTag', 'manuscriptInterval', ...
    'intervalMae', 'intervalMaeSem', 'intervalR2', 'intervalBioyield', 'intervalPremature', ...
    'intervalSafeStopRate', 'intervalRelErr', 'intervalRelErrSem', 'intervalStopMae', ...
    'intervalEp', 'intervalAlpha', 'intervalAlphaMin', 'intervalAlphaMax', ...
    'bestInterval', 'bestIntervalLabel', ...
    'bestMae', 'bestMaeSem', 'bestR2', 'bestBioyield', 'bestPremature', ...
    'bestSafeStopRate', 'bestAlpha'});
end

function m = metricsAtKey(summary, early, alphaTbl, key)
hit = summary(string(summary.krMethodKey) == string(key), :);
if isempty(hit)
    error("export_review_comparison_tables:MissingKey", "Missing interval %s.", key);
end
m = struct();
m.mae = hit.finalUpdateMae(1);
m.maeSem = col1(hit, "finalUpdateMae_sem");
m.r2 = col1(hit, "finalUpdateR2");
m.bioyield = col1(hit, "nSafeStopFail");
m.safeStopRate = col1(hit, "safeStopRate");
m.relErr = col1(hit, "relativeFinalUpdateError_mean");
m.relErrSem = col1(hit, "relativeFinalUpdateError_sem");
m.stopMae = col1(hit, "stopMae_success");
m.premature = nan;
if ~isempty(early)
    e = early(string(early.krMethodKey) == string(key), :);
    if ~isempty(e)
        m.premature = e.nEarlyStop(1);
    end
end
m.ep = nan;
m.alpha = col1(hit, "alphaDesign");
m.alphaMin = nan;
m.alphaMax = nan;
if ~isempty(alphaTbl)
    a = alphaTbl(string(alphaTbl.krMethodKey) == string(key), :);
    if ~isempty(a)
        m.ep = col1(a, "ep");
        m.alpha = col1(a, "alphaDesign");
        m.alphaMin = col1(a, "alphaDesignMin");
        m.alphaMax = col1(a, "alphaDesignMax");
    end
end
end

function T = harvestComparisonTable(cfg, gammaTag, lockedKey)
specs = { ...
    "all days (manuscript)", "jeffreys_bi_iqr15", cfg; ...
    "day 1 (2026-04-23)", "sens_harvest_0423", applyJeffreysHarvestDayPaths(cfg, "0423"); ...
    "day 2 (2026-04-29)", "sens_harvest_0429", applyJeffreysHarvestDayPaths(cfg, "0429")};
rows = cell(size(specs, 1), 1);
for i = 1:size(specs, 1)
    setting = string(specs{i, 1});
    tag = string(specs{i, 2});
    cohortCfg = specs{i, 3};
    cohort = loadStudyCohort(cohortCfg, struct("useOutlierFilter", false));
    y = double(cohort.y(:));
    k = stiffnessColumn(cohort, lockedKey, "chord");
    y = y(isfinite(y));
    kOk = k(isfinite(k));

    post = readIfExists(fullfile(cfg.out.q1, tag, "kr_methods_loocv_summary.csv"));
    post = post(string(post.methodType) == "force_abs", :);
    if ismember("variant", post.Properties.VariableNames)
        post = post(string(post.variant) == "chord", :);
    end
    [~, pix] = min(post.mae_loocv);
    postBestKey = string(post.krMethodKey(pix));
    postLocked = post(string(post.krMethodKey) == lockedKey, :);

    seq = sequentialComparisonRow(cfg, gammaTag, lockedKey, setting, i, tag, "settingValue");

    contrib = contributionBest(cfg, tag);
    rows{i} = table(setting, tag, numel(y), ...
        mean(y), std(y, 0), min(y), max(y), ...
        mean(kOk), std(kOk, 0), min(kOk), max(kOk), ...
        postBestKey, string(post.label(pix)), post.mae_loocv(pix), ...
        col1(post(pix, :), "mae_loocv_sem"), col1(post(pix, :), "r2_loocv"), ...
        lockedKey, col1(postLocked, "mae_loocv"), col1(postLocked, "mae_loocv_sem"), ...
        col1(postLocked, "r2_loocv"), ...
        seq.bestInterval, seq.bestMae, seq.bestMaeSem, seq.bestR2, ...
        seq.bestBioyield, seq.bestPremature, seq.bestSafeStopRate, seq.bestAlpha, ...
        seq.intervalMae, seq.intervalMaeSem, seq.intervalR2, ...
        seq.intervalBioyield, seq.intervalPremature, seq.intervalAlpha, ...
        contrib.caseId, contrib.deltaMae, contrib.qBH, ...
        'VariableNames', {'setting', 'analysisTag', 'n', ...
        'F_Y_mean', 'F_Y_sd', 'F_Y_min', 'F_Y_max', ...
        'k_mean', 'k_sd', 'k_min', 'k_max', ...
        'postBestInterval', 'postBestLabel', 'postBestMae', 'postBestMaeSem', 'postBestR2', ...
        'postManuscriptInterval', 'postIntervalMae', 'postIntervalMaeSem', 'postIntervalR2', ...
        'seqBestInterval', 'seqBestMae', 'seqBestMaeSem', 'seqBestR2', ...
        'seqBestBioyield', 'seqBestPremature', 'seqBestSafeStopRate', 'seqBestAlpha', ...
        'seqIntervalMae', 'seqIntervalMaeSem', 'seqIntervalR2', ...
        'seqIntervalBioyield', 'seqIntervalPremature', 'seqIntervalAlpha', ...
        'addedPredictorBest', 'addedPredictorDeltaMae', 'addedPredictorQ'});
end
T = vertcat(rows{:});
end

function out = contributionBest(cfg, tag)
out = struct("caseId", "", "deltaMae", nan, "qBH", nan);
track = dir(fullfile(cfg.out.q5, tag, "track_offline_*"));
track = track([track.isdir]);
if isempty(track)
    return;
end
path = fullfile(track(1).folder, track(1).name, "q5_model_comparison_offline.csv");
cmp = readIfExists(path);
if isempty(cmp) || ~ismember("deltaMae", cmp.Properties.VariableNames)
    return;
end
[~, ix] = min(cmp.deltaMae);
if ismember("model", cmp.Properties.VariableNames)
    out.caseId = string(cmp.model(ix));
else
    out.caseId = string(cmp.caseId(ix));
end
out.deltaMae = cmp.deltaMae(ix);
out.qBH = col1(cmp(ix, :), "qValueBH");
end

function k = stiffnessColumn(cohort, methodKey, krVariant)
k = nan(numel(cohort.y), 1);
krCol = resolveDeployKrColumn(cohort.predictorTable, methodKey, krVariant);
if ismember(krCol, cohort.predictorTable.Properties.VariableNames)
    k = double(cohort.predictorTable.(krCol));
end
end

function clearPreviousReviewTables(outDir)
oldNames = ["seq_grid.csv", "seq_best.csv", ...
    "harvest_seq_grid.csv", "harvest_seq_best.csv", ...
    "harvest_posttest_grid.csv", "harvest_posttest_best.csv", ...
    "harvest_descriptive.csv", "harvest_corr.csv", "harvest_models.csv"];
for i = 1:numel(oldNames)
    p = fullfile(outDir, oldNames(i));
    if isfile(p)
        delete(p);
    end
end
end

function T = readIfExists(path)
T = table();
if isfile(path)
    T = readtable(path, "TextType", "string");
end
end

function v = col1(tbl, name)
v = nan;
if ~isempty(tbl) && ismember(name, tbl.Properties.VariableNames)
    v = tbl.(name)(1);
end
end
