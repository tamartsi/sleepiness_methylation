############################
#Run Combat on the MESA methylation data
############################
library(mvtnorm)
library(nlme)
library(stats)
library(lumi)
library(sva)		
library(GWASTools)



####################################################################
#load needed data:
####################################################################
prepared_dat_folder <- "/data/linkage/MESA2/Projects/2018_sleepiness_methylation/Helper_data"
# phenbotypes
thepheno <- getobj(file.path(prepared_dat_folder, 'pheno_used_in_final_analysis.RData'))
thepheno <- as.data.frame(thepheno)



## Various probes that we would like to remove (sex probes, cross-reactive probes)
load(file.path(prepared_dat_folder, 'data_Probes_Remove.Rdata'))

## they were saved in a funky way
fin.probes<-unique(c(FinProbesRemv))
fin.probes <- unique(c(fin.probes, Moreprobes))

# load methylation data, identify DNAm sites to remove, and remove them. 
load(file.path(prepared_dat_folder, 'data_pre_ComBat_ESS_easy_readR.Rdata'))
toremove<-which(rownames(allDNAm) %in% fin.probes)
allDNAm<-allDNAm[-toremove,]
gc()

dim(allDNAm)
# [1] 399526    619

colnames(allDNAm) <- sapply(colnames(allDNAm), function(x) substr(x, 6,  nchar(x)))
all(colnames(allDNAm) == thepheno$idno)
# [1] TRUE

####################################################################
#Grabbing the variables
####################################################################
chip.meth<-as.factor(thepheno[,"chip_meth"])

pos.meth <-as.factor(thepheno[,"pos_meth"])
therace  <-as.factor(thepheno[,"race"])
thesite  <-as.factor(thepheno[,"site"])

mod1<-model.matrix(~pos.meth+therace+thesite)
mod1.b<-model.matrix(~therace+thesite)
cel.2.grab<-c ("bcell","tcell","nkcell","neutro") #the cell type
new.sex<- thepheno[,"gender1"]	#Sex variable
covar.contin <- thepheno[,c("age5c",cel.2.grab)] #Other continuous Variables
covar.contin <- cbind(covar.contin,new.sex)#Adding in sex and the interaction sex and Age
newcovar<-as.matrix(covar.contin)

mod2<-cbind(mod1,newcovar)
mod2.b<-cbind(mod1.b,newcovar)






############################################################
#Now running ComBat
#Running it twice to remove chip 
#and then position on chip
############################################################

normData<-ComBat(allDNAm,chip.meth,mod=mod2)
colnames(normData)<-colnames(allDNAm)

normData.2<-ComBat(normData,pos.meth,mod=mod2.b)
##################################################
#Now to save it out 
###################################################

colnames(normData.2)<- colnames(allDNAm)


normData.beta <- m2beta(normData.2)
filesave <-file.path(prepared_dat_folder , "data_DNAm_betavals_combat_all.Rdata")
save(normData.beta,file=filesave)


q("no")
