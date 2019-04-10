
### load all p-values from all analsyes, and plot histograms using ggplot. 

library(ggplot2)
library(GWASTools)
library(data.table)

MESA_res <- "MESA_EWAS_results/data_fin_sub_ESS_betavals_combat_all_includeSNP_PCs_race.Rdata"
chs_res_folder <- "CHS_EWAS_results/"
figures_folder <- "Figures/"

prep_dat_for_plot <- function(cohort, pop, pval){
	data.frame(cohort = cohort, 
						Population = pop,
						 pval = pval, stringsAsFactors = FALSE)
}


## first load MESA results: 

load(MESA_res)
for_plot <- prep_dat_for_plot("MESA", "All", results.ESS[, "Pval"])
for (i in 1:length(results.ESS.race)){
	
	temp <- prep_dat_for_plot("MESA", switch(i, "EA", "AA", "HA"), results.ESS.race[[i]][, "Pval"])
	for_plot <- rbind(for_plot, temp)
}


CHS <- read.csv("CHS_EWAS_results/CHS_Sleep_Whites_sqEDS_1_2_2019.csv"))
temp <- prep_dat_for_plot("CHS", "EA", CHS[, "P"])
for_plot <- rbind(for_plot, temp)

CHS <- read.csv("CHS_EWAS_results/CHS_Sleep_Blacks_sqEDS_1_2_2019.csv")
temp <- prep_dat_for_plot("CHS", "AA", CHS[, "P"])
for_plot <- rbind(for_plot, temp)

chs_meta <- getobj("CHS_EWAS_results/20190123_CHS_Sleep_sqEDS_meta.RData")
temp <- prep_dat_for_plot("CHS", "All", chs_meta[, "pval.Z"])
for_plot <- rbind(for_plot, temp)

## add the combined data: 

for (pop in c("EA", "AA", "All")){
	dat <- getobj(paste0("MESA_CHS_meta_results/20190110_meta_", pop, "_betavals_combat_all.RData"))		
	temp <- prep_dat_for_plot("CHS + MESA", pop, dat[, "pval.Z"])	
	for_plot <- rbind(for_plot, temp)	
}

colnames(for_plot)[1] <- "Study"


comp_lambda <- function(pvals){
	stats <- pchisq(pvals, df = 1, lower.tail = FALSE)
	median(stats)/pchisq(0.5, df = 1, lower.tail = FALSE)
}


for_plot <- data.table(for_plot)
lambdas <- for_plot[,round(comp_lambda(pval),2), by = c("Study", "Population")]
setnames(lambdas, "V1", "lambda")
lambdas <- data.frame(lambdas)
lambdas$x <- 0.66
lambdas$y <- 20000
lambdas$y[which(lambdas$Study == "CHS")] <- 21500
lambdas$y[which(lambdas$Study == "CHS + MESA")] <- 23000
lambdas$label <- paste0(lambdas$Study,  " inf=", lambdas$lambda)
#lambdas$label <- expression(paste(lambdas$Study,"=", lambdas$lambda))












p <- ggplot(for_plot,  aes(pval, fill = Study, colour = Study))
			  			

p + 
	geom_histogram( alpha = 0.4,  position = "identity") + 
	geom_text(data = lambdas, mapping = aes(x = x, y = y, label = label)) + 
	facet_wrap(~Population, ncol = 2) +
	scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
	scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) + 
	xlim(0,1) + labs(x = "p-value")
	

 ggsave(paste0(figures_folder, "20190318_all_histograms_xlim.pdf"))

