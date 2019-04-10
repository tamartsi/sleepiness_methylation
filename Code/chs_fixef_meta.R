library(qvalue)
library(data.table)
source("Code/fixef_meta_analysis.R", chdir = TRUE)


chs.res <- list(AA = "CHS_EWAS_results/CHS_Sleep_Blacks_sqEDS_1_2_2019.csv", 
				EA = "CHS_EWAS_results/CHS_Sleep_Whites_sqEDS_1_2_2019.csv")
				
chs.AA <- fread(chs.res$AA)
chs.EA <- fread(chs.res$EA)

setnames(chs.AA, "V1", "probeID")
setnames(chs.EA, "V1", "probeID")

chs <- merge(chs.AA, chs.EA, by = "probeID", suffix = c("_AA", "_EA"))
chs <- as.data.frame(chs)
rownames(chs) <- chs$probeID
meta.chs <-  metaAnalysis(betas = chs[,paste0("beta", c("_AA", "_EA"))], stdErrors= chs[,paste0("SE", c("_AA", "_EA"))], ns = chs[,paste0("N", c("_AA", "_EA"))])
meta.chs$qval <- qvalue(meta.chs$pval.Z)$qv
min(meta.chs$qval) # 0.9999985
min(meta.chs$pval.Z) # 9.455073e-06

save(meta.chs, file = "CHS_EWAS_results/20190123_CHS_Sleep_sqEDS_meta.RData")