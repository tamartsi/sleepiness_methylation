

require(GWASTools)
require(data.table)
require(ggplot2)
library(tidyr)
library(plyr)
#library(gridExtra)

load('MESA_EWAS_results/data_fin_sub_ESS_betavals_combat_all_includeSNP_PCs_race.Rdata')   ## MESA results
phen <- getobj("pheno_used_in_final_analysis.RData") ## load phenotypes

	groups <-  c("AA", "EA", "HA")
	
	list_groups = list(AA = "MESA-blacks", EA = "MESA-whites", HA = "MESA-hispanics") # this mathces the group names we use here to the ones in the MESA phenotype file. 
	
	
	meth.dat <- getobj('data_DNAm_betavals_combat_all_1.Rdata') ## not available on public repository!

	
	ids.in.analysis <- colnames(meth.dat)
	cpg.ids <- getobj('Follow_up_results/20190116_MESA_top_probes.RData')
	
	meth.dat <- meth.dat[which(rownames(meth.dat) %in% cpg.ids), ]
	meth.dat <- data.frame(meth.dat)
	meth.dat$TargetID <- rownames(meth.dat)
	
	
	list.meth.group <- vector(mode = "list", length = length(groups))
	names(list.meth.group) <- groups
	
	## for each group, create a table with methylation beta values by site by person. 
	for (i in 1:length(groups)){
		cur.group <- groups[i]
		
		## subset to phenotype of people in the group.
		group.phen <- phen[phen$race %in% list_groups[cur.group],]
		
		meth.dat.group <- meth.dat[,c("TargetID", group.phen$idno)]
		meth.l.group <- gather(meth.dat.group, "Person", "Meth_value", group.phen$idno)
		meth.l.group$Group <- cur.group
		meth.l.group$ESS <- group.phen[, "ESS"][match(meth.l.group$Person, group.phen$idno)]
		list.meth.group[[cur.group]] <- meth.l.group
		
		
	}
	meth.for.plot <- do.call(rbind, list.meth.group)
	



p <- ggplot(meth.for.plot, aes(Meth_value, ESS,  colour = Group, shape = Group)) + 
			labs(x = "Methylation (beta-values)", y  = "Sleepiness (Epworth Sleepiness Scale)")  + 
			scale_color_manual(values=c( "#E69F00", "#56B4E9",  "#999999")) +
			geom_point() + 
			scale_shape(solid = TRUE)
p + geom_smooth(method=lm, se=TRUE, fullrange=TRUE) + facet_wrap(~TargetID, nrow = 2, scales = "free") + theme_bw()
ggsave(paste0("Figures/20190116_figure2_solid.pdf"), height = 7, width = 7)

