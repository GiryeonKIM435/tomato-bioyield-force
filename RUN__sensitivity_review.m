%% =====================================================================
%  RUN__sensitivity_review
%  Reviewer sensitivity: k-update point gate, alpha quantile, harvest day.
%  Does not overwrite jeffreys_bi_iqr15 products.
%
%  Outputs: outputs/sec4_sensitivity/review/
%  =====================================================================

doKupdate = true;     % sequential replay, minBandPoints = 2 and 20 (5 is the manuscript)
doQuantile = true;    % sequential replay, p = 0.90 and 0.99 (0.95 is the manuscript)
doHarvest = true;     % refit 4.2/4.3/4.4 on each harvest day
doCollect = true;     % comparison CSVs (also splits the pooled model by day)

reuseExisting = true;

setup_paths();
cfg = PaperStudyConfig();

if doKupdate
    run_sensitivity_kupdate(cfg, struct("reuseExisting", reuseExisting));
end
if doQuantile
    run_sensitivity_quantile(cfg, struct("reuseExisting", reuseExisting));
end
if doHarvest
    run_sensitivity_harvest_review(cfg, struct("reuseExisting", reuseExisting));
end
if doCollect
    export_review_comparison_tables(cfg);
end

fprintf("RUN__sensitivity_review finished.\n");
