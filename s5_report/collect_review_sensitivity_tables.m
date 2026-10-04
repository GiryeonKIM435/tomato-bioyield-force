function results = collect_review_sensitivity_tables(cfg)
%COLLECT_REVIEW_SENSITIVITY_TABLES Comparison CSVs for the three review factors
%
% Reads manuscript and sensitivity-arm products. Does not refit models.
% Stratified harvest rows split the pooled per-sample errors by id<=50.

if nargin < 1 || isempty(cfg)
    cfg = PaperStudyConfig();
end

outDir = fullfile(cfg.out.sensitivity, "review");
if ~isfolder(outDir)
    mkdir(outDir);
end

gammaTag = "gamma_" + strrep(sprintf("%.1f", cfg.q7.gammaValues(1)), ".", "p");
primary = "jeffreys_bi_iqr15";
maxId = cfg.paper.harvestBatchAMaxId;

seqLevels = { ...
    "kupdate", "n5", primary; ...
    "kupdate", "n2", "sens_kupdate_n2"; ...
    "quantile", "p95", primary; ...
    "quantile", "p90", "sens_quantile_p90"; ...
    "quantile", "p99", "sens_quantile_p99"};

[seqGrid, seqBest] = collectSeqLevels(cfg, seqLevels, gammaTag);
writetable(seqGrid, fullfile(outDir, "seq_grid.csv"));
writetable(seqBest, fullfile(outDir, "seq_best.csv"));

harvestSeqLevels = { ...
    "harvest", "pooled", primary; ...
    "harvest", "refit_0423", "sens_harvest_0423"; ...
    "harvest", "refit_0429", "sens_harvest_0429"};
[hSeqGrid, hSeqBest] = collectSeqLevels(cfg, harvestSeqLevels, gammaTag);

postLevels = { ...
    "harvest", "pooled", primary; ...
    "harvest", "refit_0423", "sens_harvest_0423"; ...
    "harvest", "refit_0429", "sens_harvest_0429"};
[postGrid, postBest] = collectPosttestLevels(cfg, postLevels);

[hSeqGrid, hSeqBest, postGrid, postBest] = appendStratifiedHarvest( ...
    cfg, primary, gammaTag, maxId, hSeqGrid, hSeqBest, postGrid, postBest);

writetable(hSeqGrid, fullfile(outDir, "harvest_seq_grid.csv"));
writetable(hSeqBest, fullfile(outDir, "harvest_seq_best.csv"));
writetable(postGrid, fullfile(outDir, "harvest_posttest_grid.csv"));
writetable(postBest, fullfile(outDir, "harvest_posttest_best.csv"));

desc = collectHarvestDescriptive(cfg, maxId);
writetable(desc, fullfile(outDir, "harvest_descriptive.csv"));

corr = collectHarvestCorr(cfg);
writetable(corr, fullfile(outDir, "harvest_corr.csv"));

models = collectHarvestModels(cfg);
writetable(models, fullfile(outDir, "harvest_models.csv"));

results = struct();
results.outDir = outDir;
results.seqGrid = seqGrid;
results.seqBest = seqBest;
fprintf("Review sensitivity tables: %s\n", outDir);
end

function [grid, best] = collectSeqLevels(cfg, levels, gammaTag)
gridParts = cell(size(levels, 1), 1);
bestParts = cell(size(levels, 1), 1);
for i = 1:size(levels, 1)
    factor = string(levels{i, 1});
    level = string(levels{i, 2});
    tag = string(levels{i, 3});
    [gridParts{i}, bestParts{i}] = seqTablesForTag(cfg, factor, level, tag, gammaTag);
end
grid = vertcat(gridParts{:});
best = vertcat(bestParts{:});
end

function [grid, best] = seqTablesForTag(cfg, factor, level, tag, gammaTag)
summary = readIfExists(fullfile(cfg.out.q7, tag, gammaTag, "q7_design_deploy_summary.csv"));
early = readIfExists(fullfile(cfg.out.q7, tag, gammaTag, "table_online_early_stop_by_method.csv"));
pairs = readIfExists(fullfile(cfg.out.q7, tag, gammaTag, "q7_design_vs_best_pairs_force_abs.csv"));
alphaTbl = readIfExists(fullfile(cfg.out.q7, tag, gammaTag, "q7_design_alpha_by_method.csv"));
if isempty(summary)
    warning("collect_review_sensitivity_tables:MissingSeq", ...
        "Missing sequential summary for %s/%s (%s).", factor, level, tag);
    grid = emptySeqGrid();
    best = emptySeqBest();
    return;
end
summary = summary(string(summary.methodType) == "force_abs", :);
[bestKey, bestLabel] = bestKeyByMae(summary, "finalUpdateMae", "krMethodKey", "label");
bh = bhMapFromPairs(pairs, bestKey);

n = height(summary);
factorCol = repmat(factor, n, 1);
levelCol = repmat(level, n, 1);
mae = summary.finalUpdateMae;
sem = colOrNan(summary, "finalUpdateMae_sem");
r2 = colOrNan(summary, "finalUpdateR2");
bio = colOrNan(summary, "nSafeStopFail");
prem = lookupByKey(early, summary.krMethodKey, "nEarlyStop");
isMin = string(summary.krMethodKey) == bestKey;
bhFlag = false(n, 1);
for i = 1:n
    key = string(summary.krMethodKey(i));
    if key == bestKey
        bhFlag(i) = true;
    elseif isKey(bh, char(key))
        bhFlag(i) = bh(char(key));
    end
end
grid = table(factorCol, levelCol, string(summary.krMethodKey), string(summary.label), ...
    mae, sem, r2, bio, prem, isMin, bhFlag, ...
    'VariableNames', {'factor', 'level', 'krMethodKey', 'label', ...
    'seq_mae', 'seq_maeSem', 'seq_r2', 'seq_bioyield', 'seq_premature', ...
    'seq_isMin', 'seq_bhVsMin'});

best = seqBestFromSummary(factor, level, summary, early, alphaTbl, bestKey, bestLabel);
end

function best = seqBestFromSummary(factor, level, summary, early, alphaTbl, bestKey, bestLabel)
row = summary(string(summary.krMethodKey) == bestKey, :);
if isempty(row)
    best = emptySeqBest();
    return;
end
prem = lookupByKey(early, bestKey, "nEarlyStop");
if isempty(prem) || ~isfinite(prem)
    prem = nan;
else
    prem = prem(1);
end
epMed = nan; aMed = nan; aMin = nan; aMax = nan;
if ~isempty(alphaTbl)
    aRow = alphaTbl(string(alphaTbl.krMethodKey) == bestKey, :);
    if ~isempty(aRow)
        epMed = getColFirst(aRow, "ep");
        aMed = getColFirst(aRow, "alphaDesign");
        aMin = getColFirst(aRow, "alphaDesignMin");
        aMax = getColFirst(aRow, "alphaDesignMax");
    end
end
nFail = getColFirst(row, "nSafeStopFail");
nEval = getColFirst(row, "nEvaluated");
nCohort = getColFirst(row, "nCohort");
rate = getColFirst(row, "safeStopRate");
best = table(factor, level, bestKey, bestLabel, ...
    getColFirst(row, "finalUpdateMae"), getColFirst(row, "finalUpdateMae_sem"), ...
    getColFirst(row, "finalUpdateR2"), nFail, prem, rate, nEval, nCohort, ...
    getColFirst(row, "stopMae_success"), ...
    getColFirst(row, "relativeFinalUpdateError_mean"), ...
    getColFirst(row, "relativeFinalUpdateError_sem"), ...
    epMed, aMed, aMin, aMax, ...
    'VariableNames', {'factor', 'level', 'seq_bestKey', 'seq_bestLabel', ...
    'seq_bestMae', 'seq_bestMaeSem', 'seq_bestR2', 'nBioyield', 'nPremature', ...
    'safeStopRate', 'nEvaluated', 'nCohort', 'stopMae_success', ...
    'relErr_mean', 'relErr_sem', 'ep_median', 'alpha_median', 'alpha_min', 'alpha_max'});
end

function [grid, best] = collectPosttestLevels(cfg, levels)
gridParts = cell(size(levels, 1), 1);
bestParts = cell(size(levels, 1), 1);
for i = 1:size(levels, 1)
    factor = string(levels{i, 1});
    level = string(levels{i, 2});
    tag = string(levels{i, 3});
    summary = readIfExists(fullfile(cfg.out.q1, tag, "kr_methods_loocv_summary.csv"));
    pairs = readIfExists(fullfile(cfg.out.q1, tag, "kr_methods_vs_best_pairs_force_abs.csv"));
    if isempty(summary)
        warning("collect_review_sensitivity_tables:MissingPost", ...
            "Missing post-test summary for %s.", tag);
        gridParts{i} = emptyPostGrid();
        bestParts{i} = emptyPostBest();
        continue;
    end
    summary = summary(string(summary.methodType) == "force_abs", :);
    if ismember("variant", summary.Properties.VariableNames)
        summary = summary(string(summary.variant) == "chord", :);
    end
    [gridParts{i}, bestParts{i}] = postGridFromSummary(factor, level, summary, pairs);
end
grid = vertcat(gridParts{:});
best = vertcat(bestParts{:});
end

function [grid, best] = postGridFromSummary(factor, level, summary, pairs)
[bestKey, bestLabel] = bestKeyByMae(summary, "mae_loocv", "krMethodKey", "label");
bh = bhMapFromPairs(pairs, bestKey);
n = height(summary);
isMin = string(summary.krMethodKey) == bestKey;
bhFlag = false(n, 1);
for i = 1:n
    key = string(summary.krMethodKey(i));
    if key == bestKey
        bhFlag(i) = true;
    elseif isKey(bh, char(key))
        bhFlag(i) = bh(char(key));
    end
end
grid = table(repmat(factor, n, 1), repmat(level, n, 1), ...
    string(summary.krMethodKey), string(summary.label), ...
    summary.mae_loocv, colOrNan(summary, "mae_loocv_sem"), colOrNan(summary, "r2_loocv"), ...
    isMin, bhFlag, ...
    'VariableNames', {'factor', 'level', 'krMethodKey', 'label', ...
    'post_mae', 'post_maeSem', 'post_r2', 'post_isMin', 'post_bhVsMin'});
bestMae = nan; bestSem = nan; bestR2 = nan;
row = summary(string(summary.krMethodKey) == bestKey, :);
if ~isempty(row)
    bestMae = row.mae_loocv(1);
    bestSem = getColFirst(row, "mae_loocv_sem");
    bestR2 = getColFirst(row, "r2_loocv");
end
best = table(factor, level, bestKey, bestLabel, bestMae, bestSem, bestR2, ...
    'VariableNames', {'factor', 'level', 'post_bestKey', 'post_bestLabel', ...
    'post_bestMae', 'post_bestMaeSem', 'post_bestR2'});
end

function [seqGrid, seqBest, postGrid, postBest] = appendStratifiedHarvest( ...
    cfg, primary, gammaTag, maxId, seqGrid, seqBest, postGrid, postBest)

perPath = fullfile(cfg.out.q7, primary, gammaTag, "q7_design_deploy_per_sample.csv");
per = readIfExists(perPath);
sumPath = fullfile(cfg.out.q7, primary, gammaTag, "q7_design_deploy_summary.csv");
summary = readIfExists(sumPath);
if ~isempty(per) && ~isempty(summary)
    summary = summary(string(summary.methodType) == "force_abs", :);
    dayLevels = ["pooled_stratified_0423", "pooled_stratified_0429"];
    dayEarly = [true, false];
    for i = 1:2
        maskIds = stratifiedIdMask(per.id, maxId, dayEarly(i));
        sub = per(maskIds, :);
        [g, b] = seqGridFromPerSample(summary, sub, "harvest", dayLevels(i));
        seqGrid = [seqGrid; g]; %#ok<AGROW>
        seqBest = [seqBest; b]; %#ok<AGROW>
    end
else
    warning("collect_review_sensitivity_tables:NoPerSample", ...
        "Pooled sequential per-sample table missing; stratified 4.3 skipped.");
end

matPath = fullfile(cfg.out.q1, primary, "kr_benchmark_results.mat");
if isfile(matPath)
    s = load(matPath, "results");
    if isfield(s, "results") && isfield(s.results, "cvResults") && isfield(s.results, "cohort")
        dayLevels = ["pooled_stratified_0423", "pooled_stratified_0429"];
        dayEarly = [true, false];
        for i = 1:2
            [g, b] = postGridFromBenchmark(s.results, maxId, dayEarly(i), ...
                "harvest", dayLevels(i));
            postGrid = [postGrid; g]; %#ok<AGROW>
            postBest = [postBest; b]; %#ok<AGROW>
        end
    end
else
    warning("collect_review_sensitivity_tables:NoBenchmark", ...
        "Pooled LOOCV mat missing; stratified 4.2 skipped.");
end
end

function [grid, best] = seqGridFromPerSample(summary, per, factor, level)
keys = string(summary.krMethodKey);
n = numel(keys);
mae = nan(n, 1);
sem = nan(n, 1);
r2 = nan(n, 1);
bio = nan(n, 1);
prem = nan(n, 1);
labels = string(summary.label);
absByKey = cell(n, 1);
for i = 1:n
    ps = per(string(per.krMethodKey) == keys(i), :);
    if isempty(ps)
        continue;
    end
    yTrue = double(ps.yTrue);
    yHat = double(ps.y_hat_finalUpdate);
    m = calcMetrics(yTrue, yHat);
    mae(i) = m.mae;
    sem(i) = absoluteErrorSem(yTrue, yHat);
    r2(i) = m.r2;
    bio(i) = sum(string(ps.outcome) ~= "success");
    bandHigh = double(summary.gridStart(i)) + double(summary.gridWidth(i));
    fFinal = double(ps.F_finalUpdate);
    prem(i) = sum(isfinite(fFinal) & (fFinal < bandHigh - 0.5));
    absByKey{i} = abs(double(ps.finalUpdateErrorN));
end
[~, ix] = min(mae);
bestKey = keys(ix);
bhFlag = bhFlagsFromAbsErrors(absByKey, ix);
isMin = keys == bestKey;
grid = table(repmat(string(factor), n, 1), repmat(string(level), n, 1), ...
    keys, labels, mae, sem, r2, bio, prem, isMin, bhFlag, ...
    'VariableNames', {'factor', 'level', 'krMethodKey', 'label', ...
    'seq_mae', 'seq_maeSem', 'seq_r2', 'seq_bioyield', 'seq_premature', ...
    'seq_isMin', 'seq_bhVsMin'});

psBest = per(string(per.krMethodKey) == bestKey, :);
best = stratifiedSeqBest(factor, level, summary, psBest, bestKey, labels(ix), mae(ix), sem(ix), r2(ix), bio(ix), prem(ix));
end

function best = stratifiedSeqBest(factor, level, summary, ps, bestKey, bestLabel, mae, sem, r2, bio, prem)
nEval = nnz(isfinite(ps.finalUpdateErrorN));
nCohort = height(ps);
rate = mean(string(ps.outcome) == "success");
ok = string(ps.outcome) == "success" & isfinite(ps.stopErrorN);
if any(ok)
    stopMae = mean(abs(double(ps.stopErrorN(ok))));
else
    stopMae = nan;
end
rel = double(ps.relativeFinalUpdateError);
rel = rel(isfinite(rel));
if isempty(rel)
    relMean = nan;
    relSem = nan;
else
    relMean = mean(rel);
    relSem = continuousSem(rel);
end
a = double(ps.alphaDesign);
a = a(isfinite(a));
if isempty(a)
    aMed = nan; aMin = nan; aMax = nan;
else
    aMed = median(a);
    aMin = min(a);
    aMax = max(a);
end
epMed = nan;
if ismember("ep", summary.Properties.VariableNames)
    row = summary(string(summary.krMethodKey) == bestKey, :);
    if ~isempty(row)
        epMed = row.ep(1);
    end
end
best = table(string(factor), string(level), bestKey, bestLabel, mae, sem, r2, ...
    bio, prem, rate, nEval, nCohort, stopMae, relMean, relSem, ...
    epMed, aMed, aMin, aMax, ...
    'VariableNames', {'factor', 'level', 'seq_bestKey', 'seq_bestLabel', ...
    'seq_bestMae', 'seq_bestMaeSem', 'seq_bestR2', 'nBioyield', 'nPremature', ...
    'safeStopRate', 'nEvaluated', 'nCohort', 'stopMae_success', ...
    'relErr_mean', 'relErr_sem', 'ep_median', 'alpha_median', 'alpha_min', 'alpha_max'});
end

function [grid, best] = postGridFromBenchmark(benchmark, maxId, dayEarly, factor, level)
ids = double(benchmark.cohort.ids(:));
keys = string(benchmark.methodKeys(:));
variants = string(benchmark.q1Variants(:));
vi = find(variants == "chord", 1);
if isempty(vi)
    vi = 1;
end
summary = benchmark.summaryTable;
summary = summary(string(summary.methodType) == "force_abs", :);
if ismember("variant", summary.Properties.VariableNames)
    summary = summary(string(summary.variant) == "chord", :);
end
useKeys = string(summary.krMethodKey);
n = numel(useKeys);
mae = nan(n, 1);
sem = nan(n, 1);
r2 = nan(n, 1);
labels = strings(n, 1);
absByKey = cell(n, 1);
idMaskAll = stratifiedIdMask(ids, maxId, dayEarly);
for i = 1:n
    mi = find(keys == useKeys(i), 1);
    labels(i) = string(summary.label(i));
    if isempty(mi) || vi > size(benchmark.cvResults, 2)
        continue;
    end
    cv = benchmark.cvResults{mi, vi};
    if ~isstruct(cv) || ~isfield(cv, "yTrue") || numel(cv.yTrue) ~= numel(ids)
        continue;
    end
    yTrue = double(cv.yTrue(idMaskAll));
    yPred = double(cv.yPred(idMaskAll));
    m = calcMetrics(yTrue, yPred);
    mae(i) = m.mae;
    sem(i) = absoluteErrorSem(yTrue, yPred);
    r2(i) = m.r2;
    absByKey{i} = abs(yTrue - yPred);
end
[~, ix] = min(mae);
bestKey = useKeys(ix);
bhFlag = bhFlagsFromAbsErrors(absByKey, ix);
isMin = useKeys == bestKey;
grid = table(repmat(string(factor), n, 1), repmat(string(level), n, 1), ...
    useKeys, labels, mae, sem, r2, isMin, bhFlag, ...
    'VariableNames', {'factor', 'level', 'krMethodKey', 'label', ...
    'post_mae', 'post_maeSem', 'post_r2', 'post_isMin', 'post_bhVsMin'});
best = table(string(factor), string(level), bestKey, labels(ix), mae(ix), sem(ix), r2(ix), ...
    'VariableNames', {'factor', 'level', 'post_bestKey', 'post_bestLabel', ...
    'post_bestMae', 'post_bestMaeSem', 'post_bestR2'});
end

function desc = collectHarvestDescriptive(cfg, maxId)
specs = { ...
    "pooled", cfg; ...
    "refit_0423", applyJeffreysHarvestDayPaths(cfg, "0423"); ...
    "refit_0429", applyJeffreysHarvestDayPaths(cfg, "0429")};
rows = cell(0, 1);
pY = nan;
pK = nan;
for i = 1:size(specs, 1)
    level = string(specs{i, 1});
    c = specs{i, 2};
    try
        cohort = loadStudyCohort(c, struct("useOutlierFilter", false));
    catch ME
        warning("collect_review_sensitivity_tables:Cohort", "%s: %s", level, ME.message);
        continue;
    end
    y = double(cohort.y(:));
    ids = double(cohort.ids(:));
    k = stiffnessColumn(cohort, "force_s05_w30", "chord");
        rows{end + 1, 1} = descriptiveRow("harvest", level, y, k); %#ok<AGROW>
    if level == "pooled"
        m1 = ids <= maxId;
        m2 = ids > maxId;
        pY = safeRanksum(y(m1), y(m2));
        pK = safeRanksum(k(m1), k(m2));
    end
end
desc = vertcat(rows{:});
desc.p_F_Y = repmat(pY, height(desc), 1);
desc.p_k = repmat(pK, height(desc), 1);
end

function row = descriptiveRow(factor, level, y, k)
y = y(isfinite(y));
k = k(isfinite(k));
row = table(string(factor), string(level), numel(y), ...
    mean(y), std(y, 0), min(y), max(y), ...
    mean(k), std(k, 0), min(k), max(k), ...
    'VariableNames', {'factor', 'level', 'n', ...
    'F_Y_mean', 'F_Y_sd', 'F_Y_min', 'F_Y_max', ...
    'k_mean', 'k_sd', 'k_min', 'k_max'});
end

function corr = collectHarvestCorr(cfg)
levels = { ...
    "pooled", "jeffreys_bi_iqr15"; ...
    "refit_0423", "sens_harvest_0423"; ...
    "refit_0429", "sens_harvest_0429"};
parts = cell(size(levels, 1), 1);
for i = 1:size(levels, 1)
    path = fullfile(cfg.out.q5, levels{i, 2}, "predictor_correlation_vif.csv");
    T = readIfExists(path);
    if isempty(T)
        warning("collect_review_sensitivity_tables:MissingCorr", "Missing %s", path);
        continue;
    end
    n = height(T);
    rhoStars = strings(n, 1);
    partialStars = strings(n, 1);
    for r = 1:n
        rhoStars(r) = significanceStars(T.spearmanP(r));
        partialStars(r) = significanceStars(T.partialP(r));
    end
    parts{i} = table(repmat("harvest", n, 1), repmat(string(levels{i, 1}), n, 1), ...
        string(T.predictor), T.spearmanRho, T.spearmanCiLo, T.spearmanCiHi, rhoStars, ...
        T.partialR, T.partialCiLo, T.partialCiHi, partialStars, T.vif, ...
        'VariableNames', {'factor', 'level', 'predictor', ...
        'rho', 'rhoCiLo', 'rhoCiHi', 'rhoStars', ...
        'partialR', 'partialCiLo', 'partialCiHi', 'partialStars', 'vif'});
end
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    corr = table();
else
    corr = vertcat(parts{:});
end
end

function models = collectHarvestModels(cfg)
levels = { ...
    "pooled", "jeffreys_bi_iqr15"; ...
    "refit_0423", "sens_harvest_0423"; ...
    "refit_0429", "sens_harvest_0429"};
parts = cell(size(levels, 1), 1);
for i = 1:size(levels, 1)
    tag = levels{i, 2};
    track = dir(fullfile(cfg.out.q5, tag, "track_offline_*"));
    track = track([track.isdir]);
    if isempty(track)
        warning("collect_review_sensitivity_tables:MissingModels", "No 4.4 track for %s", tag);
        continue;
    end
    trackDir = fullfile(track(1).folder, track(1).name);
    summary = readIfExists(fullfile(trackDir, "offline_loocv_summary.csv"));
    cmp = readIfExists(fullfile(trackDir, "q5_model_comparison_offline.csv"));
    if isempty(summary)
        continue;
    end
    n = height(summary);
    delta = nan(n, 1);
    pW = nan(n, 1);
    q = nan(n, 1);
    if ~isempty(cmp) && ismember("model", cmp.Properties.VariableNames)
        for r = 1:n
            cid = string(summary.caseId(r));
            hit = cmp(string(cmp.model) == cid, :);
            if isempty(hit)
                continue;
            end
            delta(r) = getColFirst(hit, "deltaMae");
            pW(r) = getColFirst(hit, "pWilcoxon");
            q(r) = getColFirst(hit, "qValueBH");
        end
    end
    parts{i} = table(repmat("harvest", n, 1), repmat(string(levels{i, 1}), n, 1), ...
        string(summary.caseId), summary.r2_loocv, summary.mae_loocv, summary.n, ...
        delta, pW, q, ...
        'VariableNames', {'factor', 'level', 'caseId', 'r2', 'mae', 'n', ...
        'deltaMae', 'pWilcoxon', 'qBH'});
end
parts = parts(~cellfun(@isempty, parts));
if isempty(parts)
    models = table();
else
    models = vertcat(parts{:});
end
end

function bh = bhMapFromPairs(pairs, bestKey)
bh = containers.Map("KeyType", "char", "ValueType", "logical");
if isempty(pairs) || ~ismember("qValueBH", pairs.Properties.VariableNames)
    return;
end
sub = pairs;
if ismember("referenceMethod", sub.Properties.VariableNames)
    sub = sub(string(sub.referenceMethod) == bestKey, :);
end
for i = 1:height(sub)
    key = char(string(sub.comparisonMethod(i)));
    q = sub.qValueBH(i);
    bh(key) = isfinite(q) && q >= 0.05;
end
end

function flags = bhFlagsFromAbsErrors(absByKey, bestIdx)
n = numel(absByKey);
flags = false(n, 1);
flags(bestIdx) = true;
ref = absByKey{bestIdx};
pVals = nan(n, 1);
idx = zeros(0, 1);
rawP = zeros(0, 1);
for i = 1:n
    if i == bestIdx || isempty(absByKey{i}) || isempty(ref)
        continue;
    end
    pVals(i) = pairedAbsoluteErrorWilcoxon(absByKey{i}, ref);
    if isfinite(pVals(i))
        idx(end + 1, 1) = i; %#ok<AGROW>
        rawP(end + 1, 1) = pVals(i); %#ok<AGROW>
    end
end
if isempty(rawP)
    return;
end
fdr = applyBenjaminiHochberg(rawP);
for k = 1:numel(idx)
    flags(idx(k)) = isfinite(fdr.qValue(k)) && fdr.qValue(k) >= 0.05;
end
end

function [key, label] = bestKeyByMae(summary, maeCol, keyCol, labelCol)
vals = summary.(maeCol);
[~, ix] = min(vals);
key = string(summary.(keyCol)(ix));
label = "";
if ismember(labelCol, summary.Properties.VariableNames)
    label = string(summary.(labelCol)(ix));
end
end

function v = lookupByKey(tbl, keys, col)
keys = string(keys);
if isempty(tbl) || ~ismember(col, tbl.Properties.VariableNames) ...
        || ~ismember("krMethodKey", tbl.Properties.VariableNames)
    v = nan(numel(keys), 1);
    return;
end
v = nan(numel(keys), 1);
for i = 1:numel(keys)
    hit = tbl(string(tbl.krMethodKey) == keys(i), :);
    if ~isempty(hit)
        v(i) = hit.(col)(1);
    end
end
end

function mask = stratifiedIdMask(ids, maxId, dayEarly)
ids = double(ids(:));
if dayEarly
    mask = ids <= maxId;
else
    mask = ids > maxId;
end
end

function k = stiffnessColumn(cohort, methodKey, krVariant)
k = nan(numel(cohort.y), 1);
krCol = resolveDeployKrColumn(cohort.predictorTable, methodKey, krVariant);
if ismember(krCol, cohort.predictorTable.Properties.VariableNames)
    k = double(cohort.predictorTable.(krCol));
end
end

function p = safeRanksum(a, b)
a = a(isfinite(a));
b = b(isfinite(b));
p = nan;
if numel(a) >= 2 && numel(b) >= 2
    try
        p = ranksum(a, b);
    catch
        p = nan;
    end
end
end

function stars = significanceStars(pVal)
if ~isfinite(pVal)
    stars = "";
elseif pVal < 0.001
    stars = "***";
elseif pVal < 0.01
    stars = "**";
elseif pVal < 0.05
    stars = "*";
else
    stars = "";
end
end

function T = readIfExists(path)
T = table();
if isfile(path)
    T = readtable(path, "TextType", "string");
end
end

function v = colOrNan(tbl, name)
if ismember(name, tbl.Properties.VariableNames)
    v = tbl.(name);
else
    v = nan(height(tbl), 1);
end
end

function v = getColFirst(tbl, name)
v = nan;
if ismember(name, tbl.Properties.VariableNames) && height(tbl) > 0
    v = tbl.(name)(1);
end
end

function T = emptySeqGrid()
T = table(strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    false(0, 1), false(0, 1), ...
    'VariableNames', {'factor', 'level', 'krMethodKey', 'label', ...
    'seq_mae', 'seq_maeSem', 'seq_r2', 'seq_bioyield', 'seq_premature', ...
    'seq_isMin', 'seq_bhVsMin'});
end

function T = emptySeqBest()
T = table(strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'factor', 'level', 'seq_bestKey', 'seq_bestLabel', ...
    'seq_bestMae', 'seq_bestMaeSem', 'seq_bestR2', 'nBioyield', 'nPremature', ...
    'safeStopRate', 'nEvaluated', 'nCohort', 'stopMae_success', ...
    'relErr_mean', 'relErr_sem', 'ep_median', 'alpha_median', 'alpha_min', 'alpha_max'});
end

function T = emptyPostGrid()
T = table(strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), false(0, 1), false(0, 1), ...
    'VariableNames', {'factor', 'level', 'krMethodKey', 'label', ...
    'post_mae', 'post_maeSem', 'post_r2', 'post_isMin', 'post_bhVsMin'});
end

function T = emptyPostBest()
T = table(strings(0, 1), strings(0, 1), strings(0, 1), strings(0, 1), ...
    zeros(0, 1), zeros(0, 1), zeros(0, 1), ...
    'VariableNames', {'factor', 'level', 'post_bestKey', 'post_bestLabel', ...
    'post_bestMae', 'post_bestMaeSem', 'post_bestR2'});
end
