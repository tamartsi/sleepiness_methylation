require(GWASTools)

cpg.ids.mesa <- getobj("Follow_up_results/20190116_MESA_top_probes.RData")
cpg.ids.meta <- getobj("Follow_up_results/20190116_META_AA_top_probes.RData")

cpg.ids <- c(cpg.ids.mesa, cpg.ids.meta)


mesa.results.file <- "MESA_EWAS_results/data_fin_sub_ESS_betavals_combat_all_includeSNP_PCs_race.Rdata"

chs.EA.results.file <- "CHS_EWAS_results/CHS_Sleep_Whites_sqEDS_1_2_2019.csv"
chs.AA.results.file <- "CHS_EWAS_results/CHS_Sleep_Blacks_sqEDS_1_2_2019.csv"

chs.meta.file <- "CHS_EWAS_results/20190123_CHS_Sleep_sqEDS_meta.RData"


load(mesa.results.file)

temp <- results.ESS.race[[1]]
temp <- temp[match(cpg.ids, rownames(temp)),]
all(cpg.ids == rownames(temp))
res <- data.frame(cpg.ids, MESA_EA_beta = round(temp[,1], 2), MESA_EA_pval = formatC(temp[,5], format = "E", digits = 2))

temp <- results.ESS.race[[2]]
temp <- temp[match(cpg.ids, rownames(temp)),]
all(cpg.ids == rownames(temp))
res$MESA_AA_beta <- round(temp[,1], 2)
res$MESA_AA_pval = formatC(temp[,5], format = "E", digits = 2)

temp <- results.ESS.race[[3]]
temp <- temp[match(cpg.ids, rownames(temp)),]
all(cpg.ids == rownames(temp))
res$MESA_HA_beta <- round(temp[,1], 2)
res$MESA_HA_pval = formatC(temp[,5], format = "E", digits = 2)


temp <- results.ESS
temp <- temp[match(cpg.ids, rownames(temp)),]
all(cpg.ids == rownames(temp))
res$MESA_All_beta <- round(temp[,1], 2)
res$MESA_All_pval = formatC(temp[,5], format = "E", digits = 2)



temp <- read.csv(chs.EA.results.file, header = TRUE)
temp <- temp[match(cpg.ids, temp$X),]
all(cpg.ids == temp$X)
res$CHS_EA_beta <- round(temp$beta, 2)
res$CHS_EA_pval = formatC(temp$P, format = "E", digits = 2)


temp <- read.csv(chs.AA.results.file, header = TRUE)
temp <- temp[match(cpg.ids, temp$X),]
all(cpg.ids == temp$X)
res$CHS_AA_beta <- round(temp$beta, 2)
res$CHS_AA_pval = formatC(temp$P, format = "E", digits = 2)



temp <- getobj(chs.meta.file)
temp <- temp[match(cpg.ids, temp$probeID),]
all(cpg.ids == temp$probeID)
res$CHS_All_beta <- round(temp$beta.meta, 2)
res$CHS_All_pval = formatC(temp$pval.Z, format = "E", digits = 2)


write.csv(res, file = "Follow_up_results/20190211_supp_table2_across_groups.csv")

