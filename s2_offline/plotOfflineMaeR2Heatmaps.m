function outPaths = plotOfflineMaeR2Heatmaps(summaryTable, pairByType, outDir, cfg, maeClim)
%plotOfflineMaeR2Heatmaps 結果4.2: LOOCV MAE ヒートマップ（MAPE と R^2 を併記）
%
% 色 = LOOCV MAE、セル注記 = mean±SEM、MAPE、R^2 の順。
% MIN 赤枠 + BH 非有意（q>=0.05）ピンク枠。
% pairByType: methodType ごとのペア表 struct（force_abs / force_trailing）
%   または旧形式の単一 table（後方互換）。
% maeClim を渡せば色尺度を強制できる。

if nargin < 5
    maeClim = [];
end

if ~isfolder(outDir)
    mkdir(outDir);
end

methodTypes = "force_abs";
prefixByType = struct("force_abs", "force_abs", "force_trailing", "force_trail");
if isempty(maeClim) || numel(maeClim) ~= 2 || ~all(isfinite(maeClim))
    maeClim = computeGlobalMaeHeatmapClim(summaryTable, [], cfg);
end

outPaths = strings(0, 1);
for ti = 1:numel(methodTypes)
    mt = methodTypes(ti);
    prefix = prefixByType.(char(mt));
    sub = summaryTable(string(summaryTable.methodType) == mt, :);
    if isempty(sub)
        continue;
    end

    pairTable = resolvePairTableForType(pairByType, mt);
    refInfo = resolveKrHeatmapReference(pairTable, sub, "mae_loocv");

    [maeMat, semMat, starts, widths, ~] = buildKrGridMatrices( ...
        sub, "mae_loocv", "mae_loocv_sem", mt);
    if isempty(maeMat) || all(isnan(maeMat(:)))
        continue;
    end
    [r2Mat, ~, ~, ~, ~] = buildKrGridMatrices(sub, "r2_loocv", "", mt);
    [mapeMat, ~, ~, ~, ~] = buildKrGridMatrices(sub, "relativeError_loocv", "", mt);
    annotationLines = buildPosttestAnnotationLines(mapeMat, r2Mat, size(maeMat));

    nondiffMask = [];
    if strlength(string(refInfo.methodKey)) > 0
        nondiffMask = buildKrWilcoxonNondiffMask( ...
            pairTable, mt, refInfo.methodKey, ...
            struct("referenceVariant", refInfo.variant));
    end

    layout = resolveKrHeatmapLayout(prefix, starts, widths);
    outPath = fullfile(outDir, "fig4_2_offline_mae_r2_" + prefix + ".png");
    plotKrGridHeatmap(maeMat, semMat, starts, widths, ...
        title=sprintf("Offline LOOCV MAE (%s, %s)", ...
        strrep(prefix, "_", "\_"), string(cfg.deploy.krVariant)), ...
        outPath=outPath, cfg=cfg, clim=maeClim, colorbarLabel="LOOCV MAE [N]", ...
        highlightMode="minNondiffFdr", scaleMode="abs", showSem=true, ...
        figureSize=layout.figSize, compactText=layout.compactText, ...
        xLabel=layout.xLabel, yLabel=layout.yLabel, ...
        nondiffMask=nondiffMask, ...
        annotationLines=annotationLines);
    outPaths(end + 1, 1) = outPath; %#ok<AGROW>
end

end

function pairTable = resolvePairTableForType(pairByType, methodType)
pairTable = [];
if isempty(pairByType)
    return;
end
if istable(pairByType)
    pairTable = pairByType;
    return;
end
if ~isstruct(pairByType)
    return;
end
key = char(methodType);
if isfield(pairByType, key)
    pairTable = pairByType.(key);
elseif methodType == "force_trailing" && isfield(pairByType, "force_trail")
    pairTable = pairByType.force_trail;
end
end

function annotationLines = buildPosttestAnnotationLines(mapeMat, r2Mat, matSize)
annotationLines = cell(matSize);
annotationLines(:) = {""};
hasMape = ~isempty(mapeMat) && isequal(size(mapeMat), matSize);
hasR2 = ~isempty(r2Mat) && isequal(size(r2Mat), matSize);
for yi = 1:matSize(1)
    for xi = 1:matSize(2)
        parts = strings(0, 1);
        if hasMape && isfinite(mapeMat(yi, xi))
            parts(end + 1) = sprintf('MAPE=%.1f%%', 100 * mapeMat(yi, xi)); %#ok<AGROW>
        end
        if hasR2 && isfinite(r2Mat(yi, xi))
            parts(end + 1) = string(formatHeatmapR2Annotation(r2Mat(yi, xi))); %#ok<AGROW>
        end
        if ~isempty(parts)
            annotationLines{yi, xi} = char(strjoin(parts, newline));
        end
    end
end
end
