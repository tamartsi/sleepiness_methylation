library(qvalue)
library(data.table)
source('Code/fixef_meta_analysis.R', chdir = TRUE)


## load MESA results
load("MESA_EWAS_results/data_fin_sub_ESS_betavals_combat_all_includeSNP_PCs_race.Rdata")

chs.res <- list(AA = "CHS_EWAS_results/CHS_Sleep_Blacks_sqEDS_1_2_2019.csv", 
				EA = "CHS_EWAS_results/CHS_Sleep_Whites_sqEDS_1_2_2019.csv")
				

mesa.all <-results.ESS
colnames(mesa.all)<-c("Effect","SE","N","DF lost","Pvalue")


##################################
#European Americans results
##################################
mesa.EA <-results.ESS.race[[1]]
colnames(mesa.EA)<-c("Effect","SE","N","DF lost","Pvalue")

##################################
#African Americans results
##################################
mesa.AA <-results.ESS.race[[2]]
colnames(mesa.AA) <-c("Effect","SE","N","DF lost","Pvalue")


##################################
#Hispanic Americans results
##################################
mesa.HA <-results.ESS.race[[3]]
colnames(mesa.HA)<-c("Effect","SE","N","DF lost","Pvalue")


chs.AA <- fread(chs.res$AA)
chs.AA <- chs.AA[V1 %in% rownames(mesa.AA)]


chs.EA <- fread(chs.res$EA)
chs.EA <- chs.EA[V1 %in% rownames(mesa.EA)]

### CHS has a few less probes than MESA, so update MESA: 
all(rownames(mesa.all[match(chs.AA$V1, rownames(mesa.all)),])  == chs.AA$V1) # TRUE
all(chs.AA$V1 == chs.EA$V1) # TRUE
## therefore: 
mesa.all <- mesa.all[match(chs.AA$V1, rownames(mesa.all)),]
mesa.EA <- mesa.EA[match(chs.AA$V1, rownames(mesa.EA)),]
mesa.AA <- mesa.AA[match(chs.AA$V1, rownames(mesa.AA)),]
mesa.HA <- mesa.HA[match(chs.AA$V1, rownames(mesa.HA)),]


meta.AA <- metaAnalysis(betas = cbind(mesa.AA[,1], chs.AA$beta), stdErrors= cbind(mesa.AA[,2], chs.AA$SE), ns = cbind(mesa.AA[,3], chs.AA$N))
meta.AA$qval <- qvalue(meta.AA$pval.Z)$qv
sum(meta.AA$qval<=0.05) # 1  # 14 when considering all methylation sites. 


mesa.AA.qval <- qvalue(mesa.AA[,5])$qv

meta.EA <- metaAnalysis(betas = cbind(mesa.EA[,1], chs.EA$beta), stdErrors= cbind(mesa.EA[,2], chs.EA$SE), ns = cbind(mesa.EA[,3], chs.EA$N))
meta.EA$qval <- qvalue(meta.EA$pval.Z)$qv
min(meta.EA$qval) # 0.6621006

meta.all <-  metaAnalysis(betas = cbind(mesa.all[,1], chs.EA$beta, chs.AA$beta), stdErrors= cbind(mesa.all[,2], chs.EA$SE, chs.AA$SE), ns = cbind(mesa.all[,3], chs.EA$N, chs.AA$N))
meta.all$qval <- qvalue(meta.all$pval.Z)$qv
min(meta.all$qval) # 0.9644338



### save results from AA -- add CHS and MESA results: 
meta.AA$beta.mesa <- mesa.AA[,1]
meta.AA$se.mesa <- mesa.AA[,2]
meta.AA$pval.mesa <- mesa.AA[,5]

meta.AA$beta.chs <- chs.AA$beta
meta.AA$se.chs <- chs.AA$SE
meta.AA$pval.chs <- chs.AA$P


save(meta.EA, file = paste0("MESA_CHS_meta_results/20190110_meta_EA_betavals_combat_all.RData"))

save(meta.AA, file = paste0("MESA_CHS_meta_results/20190110_meta_AA_betavals_combat_all.RData"))

save(meta.all, file = paste0("MESA_CHS_meta_results/20190110_meta_All_betavals_combat_all.RData"))


