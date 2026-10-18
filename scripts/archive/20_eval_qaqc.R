# suppressPackageStartupMessages({ library(data.table); library(fs); library(qs) })
# source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config.R")
# 
# met_files <- dir_ls(MET_DIR, glob="metrics_*.csv", type="file", fail=FALSE)
# .log("Metric files: %d", length(met_files))
# 
# if (length(met_files)) {
#   all <- rbindlist(lapply(met_files, fread), use.names=TRUE, fill=TRUE)
#   print(all[, .(
#     n = .N,
#     mean_f1 = mean(f1, na.rm=TRUE),
#     mean_auc = mean(auc, na.rm=TRUE),
#     mean_precision = mean(precision, na.rm=TRUE),
#     mean_recall = mean(recall, na.rm=TRUE)
#   ), by=.(var_group, event)][order(var_group, event)])
# }

# 20260224 ----
suppressPackageStartupMessages({ library(data.table); library(fs); library(qs) })
source("V:/users/hquintal/phd1_cluster_southeast_v2/scripts/14_eval_config_county.R")

met_files <- dir_ls(MET_DIR_COUNTY, glob="metrics_*.csv", type="file", fail=FALSE)
.log("County metric files: %d", length(met_files))

if (length(met_files)) {
  all <- rbindlist(lapply(met_files, fread), use.names=TRUE, fill=TRUE)
  
  # Keep only ALL_COUNTIES row for quick global summary (optional)
  all_all <- all[GEOID == "ALL_COUNTIES"]
  
  print(all_all[, .(
    n = .N,
    mean_f1 = mean(f1, na.rm=TRUE),
    mean_auc = mean(auc, na.rm=TRUE),
    mean_precision = mean(precision, na.rm=TRUE),
    mean_recall = mean(recall, na.rm=TRUE)
  ), by=.(var_group, event)][order(var_group, event)])
}