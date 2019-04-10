require(GWASTools)
require(data.table)
require(ggplot2)
require(plyr)

figures_folder <- "Figures/"

fileToLoad <- 'data_DNAm_betavals_combat_all_1.Rdata' ## not available on public repository!
load(fileToLoad)



thepheno<- getobj(file.path(indata_folder, 'pheno_used_in_final_analysis.RData'))  ## not available on public repository!


EDS_IDs <- as.character(thepheno$MESA_IID[which(thepheno$ESS > 10)])   ## EDS > 10

noEDS_IDs <- as.character(thepheno$MESA_IID[which(thepheno$ESS <= 10)]) 

## bottom quartile
low_ESS_IDs <- as.character(thepheno$MESA_IID[order(thepheno$ESS)[1:(floor(nrow(thepheno)/4))]]) 

### check what the value of this bottom quartile is: 
max(thepheno$ESS[match(low_ESS_IDs, thepheno$MESA_IID)])
# [1] 3

# summary(thepheno$ESS)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 0.000   3.000   5.000   6.048   8.000  23.000


## compute mean methylation values across people:

meth_means <- colMeans(normData.beta)

for_plot <- data.frame(IDs = EDS_IDs, Group = "EDS", Meth_means = meth_means[match(EDS_IDs, names(meth_means))])
for_plot <- rbind(for_plot, 
				data.frame(IDs = noEDS_IDs, Group = "No EDS", Meth_means = meth_means[match(noEDS_IDs, names(meth_means))]), 
				data.frame(IDs = low_ESS_IDs, Group = "Low ESS", Meth_means = meth_means[match(low_ESS_IDs, names(meth_means))]))
				
				
### density plots of means: 
p <- ggplot(for_plot, aes(x = Meth_means, fill = Group)) + geom_density(alpha = 0.4) + labs(x = "Methylation means")

 ggsave(paste0(figures_folder, "20190318_mean_methylation_densities.pdf"))
 
 medians <- ddply(for_plot, "Group", summarise, Median = median(Meth_means))
 
 ### density plots of means: 
p <- ggplot(for_plot, aes(x = Meth_means, fill = Group)) + 
				geom_density(alpha = 0.4) + 
				labs(x = "Methylation means") + 
				geom_vline( data = medians, aes(xintercept=Median, colour = Group))
 ggsave(paste0(figures_folder, "20190320_mean_methylation_densities_with_median.pdf"))
				
for_plot <- data.table(for_plot)
for_plot[,mean(Meth_means), by = "Group"]


     Group        V1
1:     EDS 0.5113739
2:  No EDS 0.5113567
3: Low ESS 0.5113823


for_plot[,median(Meth_means), by = "Group"]
     Group        V1
1:     EDS 0.5114689
2:  No EDS 0.5114729
3: Low ESS 0.5115241