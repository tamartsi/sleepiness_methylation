require(GWASTools)


phen.file <- "pheno_used_in_final_analysis.RData" ## not available on the public repository!
covariate.reg.file <- 'Follow_up_results/Sensitivity/20190210_covariates_associations.csv'
cpg.ids.mesa <- getobj('Follow_up_results/20190116_MESA_top_probes.RData')
cpg.ids.meta <- getobj('Follow_up_results/20190116_META_AA_top_probes.RData')

cpg.ids <- c(cpg.ids.mesa, cpg.ids.meta)



fileToLoad <- 'data_DNAm_betavals_combat_all_1.Rdata' ## not available on the public repository!
normData.2 <- getobj(fileToLoad)
Tokeep<-which(rownames(normData.2) %in% cpg.ids)
normData.2<-normData.2[Tokeep,]



##########################################################
#Getting the phenotype data 
##########################################################
#
pheno<- getobj(phen.file)

#####################################################################################
#The residual cell type enrichment varianbles we will need
#####################################################################################
cell.vars <- c("bcell","tcell","nkcell","neutro")


#####################################################################################
#Getting the various sensitivity parameter that we want

#####################################################################################

covar.res <- read.csv(covariate.reg.file, header = TRUE, stringsAsFactors = FALSE)
sensitivity.vars <- covar.res$phenotype[which(covar.res$pval < 0.1)]
# [1] sleep_duration       whi_insomnia         AHI
# [4] per90                sleep_dur_below_5hrs CESD
# [7] total_fat            total_carbs

### add to phenotype with that name, for clarity:

pheno$sleep_duration <- pheno$avgmainsleep5
pheno$whi_insomnia <- pheno$whiirs5c
pheno$AHI <- pheno$rdi3p5
pheno$per90 <- pheno$pctlt905
pheno$sleep_dur_below_5hrs <-  ifelse(pheno$avgmainsleep5< 60*5,1,0)
pheno$CESD <- pheno$cesd5c
pheno$total_fat <- pheno$tfatn5c
pheno$total_carbs <- pheno$tcarbn5c




medication.vars <- c("apsy5c", "ntca5c", "tca5c", "tcap5c")

PCs <- paste0("PC", 1:5)


mean.impute.vec <- function(x){
	inds.mis <- which(is.na(x))
	if (length(inds.mis) > 0) x[inds.mis] <- mean(x, na.rm = TRUE)
	return(x)
}
mean.impute.all.cols <- function(mat){
	
	for (i in 1:ncol(mat)){
		mat[,i] <- mean.impute.vec(mat[,i])
	}
	return(mat)
}


### prepare results object for sensitivity analysis:
res.colnames <- c("probeID", "m1_beta", "m1_pval", paste(rep(paste0("m", 2:(length(sensitivity.vars ) + 2)), each =3), rep(c("beta", "pval", "percent"), times =5), sep = "_"), c("nomed_beta", "nomed_pval", "nomed_percent"))
res <- matrix(NA, ncol = length(res.colnames), nrow = length(cpg.ids))
res <- as.data.frame(res)
colnames(res) <- res.colnames
res$probeID <- cpg.ids

mod.n <- rep(NA, length(sensitivity.vars ) + 3)
names(mod.n) <- c(paste0("mod", 1:(length(sensitivity.vars) + 2)), "nomed")


pheno$MESA_IID <- paste0("MESA_", pheno$idno)


#### perform analysis using all people, and using AAs specifically. 
for (group in c("all", "AA")){
	cur.res <- res
	cur.n <- mod.n
	
	if (group == "all") {
		pheno.use <- pheno
		SNPPcs <- read.csv("PCs_from_Yongmei/MESA_PCA_Combined.csv")
		SNPPcs[,PCs] <- mean.impute.all.cols(SNPPcs[,PCs])
		pheno.use <- merge(pheno.use, SNPPcs[,c("idno",PCs)], by = "idno")
		base.formula <- paste("sqrt(ESS)~meth+as.factor(race)+ as.factor(gender1) +", 
								"age5c + bcell + tcell + nkcell + neutro +",
								 "PC1 + PC2 + PC3 + PC4 + PC5")
		
	}
	if (group == "AA") {
		pheno.use <- pheno[grep("blacks", pheno$race),]
		SNPPcs <- read.csv("PCs_from_Yongmei/MESA_PCA_AA.csv")
		SNPPcs[,PCs] <- mean.impute.all.cols(SNPPcs[,PCs])
		pheno.use <- merge(pheno.use, SNPPcs[,c("idno",PCs)], by = "idno")
		base.formula <- paste("sqrt(ESS)~ meth+ as.factor(gender1) +", 
								"age5c + bcell + tcell + nkcell + neutro +",
								 "PC1 + PC2 + PC3 + PC4 + PC5")
	}

	
	for (i in 1:length(cpg.ids)){
		probe <- cpg.ids[i]
		meth <- normData.2[match(probe, rownames(normData.2)),]
		meth.dat <- data.frame(meth, MESA_IID = names(meth))
		pheno.i <- merge(pheno.use, meth.dat, by = "MESA_IID")
		
		mod1 <-  lm(base.formula, data = pheno.i)
		cur.res[i,c("m1_beta", "m1_pval")] <- summary(mod1)$coef[2,c(1,4)]
		
		if (i == 1){
		cur.n["mod1"] <- nrow(model.matrix(mod1))
		}
		

		for (sens.i in 2:(length(sensitivity.vars) + 1)){
			
			mod <-  lm(paste(base.formula, 
							"+", sensitivity.vars[sens.i-1]), data = pheno.i)
			cur.res[i,c(paste0("m", sens.i, "_beta"), 
						paste0("m", sens.i, "_pval"), 
						paste0("m", sens.i, "_percent"))] <- c(summary(mod)$coef[2,c(1,4)], 
	#															summary(mod)$coef[2,1]/summary(mod1)$coef[2,1]*100)
																abs((summary(mod)$coef[2,1]-summary(mod1)$coef[2,1])/summary(mod1)$coef[2,1])*100)
																
																		
			if (i == 1){
				cur.n[paste0("mod", sens.i)] <- nrow(model.matrix(mod))
			}
		
		}
		
		
		### now a model with all sensitivity variables: 
		mod <-  lm(paste(base.formula, 
							"+" , paste0(sensitivity.vars, collapse = " + ")), data = pheno.i)
		cur.res[i,c("m10_beta", "m10_pval", "m10_percent")] <- c(summary(mod)$coef[2,c(1,4)], 
																abs((summary(mod)$coef[2,1]-summary(mod1)$coef[2,1])/summary(mod1)$coef[2,1])*100)
		
		if (i == 1){
				cur.n[paste0("mod", length(sensitivity.vars) + 2)] <- nrow(model.matrix(mod))
			}
		


			
		ind.meds <- which(rowSums(pheno.i[,medication.vars ])>1)
		mod <-  lm(base.formula, data = pheno.i[-ind.meds,])
		cur.res[i,c("nomed_beta", "nomed_pval", "nomed_percent")] <- c(summary(mod)$coef[2,c(1,4)], 
																abs((summary(mod)$coef[2,1]-summary(mod1)$coef[2,1])/summary(mod1)$coef[2,1])*100)
		
		if (i == 1){
			cur.n["nomed"] <- nrow(model.matrix(mod))
		}


		
	}
	
	write.csv(cur.res, file = paste0("Follow_up_results/Sensitivity/20190213_sensitivity_results_", group, ".csv"))	
	write.csv(cur.n, file = paste0("Follow_up_results/Sensitivity/20190213_sensitivity_samplesizes_", group, ".csv"))	
	  
}
	
	



