
library(qvalue)
require(GWASTools)

indata_folder <- "/data/linkage/MESA2/Projects/2018_sleepiness_methylation/Helper_data"
results_folder <- "/data/linkage/MESA2/Projects/2018_sleepiness_methylation/Results"
PCs_folder <- "/data/linkage/MESA2/Genotypes/PCs_from_Yongmei"

### read methylation data an dphenotypes. 

fileToLoad <- file.path(indata_folder,"data_DNAm_betavals_combat_all.Rdata" )
normData.2 <- getobj(fileToLoad)



thepheno<-getobj(file.path(indata_folder, 'pheno_used_in_final_analysis.RData'))
thepheno <- as.data.frame(thepheno)

AllPop<-thepheno$race

#####################################################################################
#Make an item to hold PC
#####################################################################################

PC.race<-vector(mode="list",length=3)
PC.all<-NULL

#####################################################################################
#Get the SNP AA PCs
#From Yongmei's group. Have to impute the ones that are missing
#Missing are four AA participants
#####################################################################################
#These are the AA PCs
SNPPcs<-read.csv(file.path(PCs_folder,"MESA_PCA_AA.csv"))

AA.SNP<-subset(SNPPcs,idno %in% thepheno$idno)
AA.pheno<-subset(thepheno,race=="MESA-blacks")

ToAdd<-AA.pheno$idno[-which(AA.pheno$idno %in% AA.SNP$idno)]
newSNP<-matrix(NA,length(ToAdd),ncol(AA.SNP)-1)


FinAdd<-data.frame(idno=ToAdd,newSNP)
colnames(FinAdd)[-1]<-colnames(AA.SNP)[-1]
PutTog<-rbind(AA.SNP,FinAdd)

PutTog<-PutTog[match(AA.pheno$idno,PutTog$idno),]

all.equal(as.character(PutTog$idno),as.character(AA.pheno$idno))



for(i in 2:51){
  PutTog[is.na(PutTog[,i]),i]<-mean(PutTog[,i],na.rm=T)
}


PC.race[[2]]<-PutTog[,-1]
rownames(PC.race[[2]])<-PutTog[,1]

#####################################################################################
#Get the SNP HA PCs
#####################################################################################
#These are the HA PCs
SNPPcs<-read.csv(file.path(PCs_folder,"MESA_PCA_HA.csv"))


HA.SNP<-subset(SNPPcs,idno %in% thepheno$idno)
HA.pheno<-subset(thepheno,race=="MESA-hispanics")


PutTog<-rbind(HA.SNP)

PutTog<-PutTog[match(HA.pheno$idno,PutTog$idno),]

all.equal(as.character(PutTog$idno),as.character(HA.pheno$idno))

PC.race[[3]]<-PutTog[,-1]
rownames(PC.race[[3]])<-PutTog[,1]
#####################################################################################
#Get the SNP EA PCs
#####################################################################################
#These are the EA PCs
SNPPcs<-read.csv(file.path(PCs_folder,"MESA_PCA_EA.csv"))

EA.SNP<-subset(SNPPcs,idno %in% thepheno$idno)
EA.pheno<-subset(thepheno,race=="MESA-whites")


PutTog<-rbind(EA.SNP)

PutTog<-PutTog[match(EA.pheno$idno,PutTog$idno),]

all.equal(as.character(PutTog$idno),as.character(EA.pheno$idno))

PC.race[[1]]<-PutTog[,-1]
rownames(PC.race[[1]])<-PutTog[,1]


#####################################################################################
#Get the SNP Combined PCs
#####################################################################################

#These are the Combined PCs
SNPPcs<-read.csv(file.path(PCs_folder,"MESA_PCA_Combined.csv"))


All.SNP<-subset(SNPPcs,idno %in% thepheno$idno)


ToAdd<-thepheno$idno[-which(thepheno$idno %in% All.SNP$idno)]

newSNP<-matrix(NA,length(ToAdd),ncol(All.SNP)-1)

FinAdd<-data.frame(idno=ToAdd,newSNP)
colnames(FinAdd)[-1]<-colnames(All.SNP)[-1]

PutTog<-rbind(All.SNP,FinAdd)



PutTog<-PutTog[match(thepheno$idno,PutTog$idno),]

all.equal(as.character(PutTog$idno),as.character(thepheno$idno))

##Fill in four missing

temp.AASNP<-subset(PutTog,thepheno$race=="MESA-blacks")
TheMeans<-colMeans(temp.AASNP[,-1],na.rm=T)

ToFill<-which(is.na(rowSums(PutTog[,-1])))
for(i in ToFill){
  PutTog[i,-1]<-TheMeans
}

PC.all<-PutTog[,-1]
rownames(PC.all)<-PutTog[,1]
#####################################################################################
#Define the race variables
#####################################################################################


zerace<-c("MESA-whites","MESA-blacks","MESA-hispanics")
names(PC.race)<-zerace


#####################################################################################
#The residual cell type enrichment varianbles we will need
#####################################################################################


cel.2.grab<- c("bcell","tcell","nkcell","neutro") #the cell type



#####################################################################################
#The function to grab the covariates that may be of interest
#####################################################################################



func.grab.cov<-function(whatpheno){
  new.sex<- whatpheno[,"gender1"]
  covar.contin<-whatpheno[,c("age5c",cel.2.grab)]
  covar.contin<-cbind(covar.contin,new.sex)
  return(as.matrix(covar.contin))
}


##########################################################################################
#Function below does fast linear regression using
#The residuals and what have you
##########################################################################################
fastReg.func<-function(Y,XX,CpG){
  numExplan<-ncol(XX)
  XXproj<-solve(t(XX)%*%XX)%*%t(XX)
  residY<-Y-XX%*%XXproj%*%Y
  onDNA<-XXproj%*%CpG
  residDNA<-CpG-XX%*%onDNA
  par1<-colSums(residDNA^2)
  par2<-colSums(residDNA*c(residY))
  Est2<-par2/par1
  PredY<-residDNA*matrix(Est2,nrow=nrow(residDNA),ncol=ncol(residDNA),byrow=T)
  newResid<-PredY-c(residY)
  SSQ<-colSums(newResid^2)
  SIGMA<-SSQ/(length(Y)-numExplan-1)
  Tstat<-Est2*sqrt(par1/SIGMA)
  Pval<-2*pt(abs(Tstat),lower.tail=F,df=length(Y)-numExplan-1)
  return(cbind(c(Est2),c(Est2/Tstat),length(Y),numExplan,Pval))
}

##########################################################################################
#This function will perform the analysis for a transformed variable Y (given a funcTransform)
#For ESS, pass in the sqrt function. 
#Run the analyses
##########################################################################################
myfunc<-function(Y,funcTransform,race=NULL){
  if(!is.null(race)){
    newrace<-which(thepheno$race==race)
    Overall.PCs<-as.matrix(PC.race[[race]])
    toKEEP<-newrace
  }
  if(is.null(race)){
    Overall.PCs<-as.matrix(PC.all)
    toKEEP<-which(!is.na(Y))
  }
  AllDNAm<-t(normData.2[,toKEEP])
  phenoKeep<-thepheno[toKEEP,]  
  therace  <-as.factor(phenoKeep[,"race"])
  site <-as.character(phenoKeep$site)
  newcovar<-func.grab.cov(phenoKeep)
  newY<-funcTransform(Y[toKEEP])
  #######################################################
  #Perform The results
  #######################################################
  N<-length(newY)
  if(!is.null(race)){
    XX1<-model.matrix(~newcovar+site+Overall.PCs[,1:5])  
  }
  if(is.null(race)){
    XX1<-model.matrix(~newcovar+site+Overall.PCs[,1:5]+therace)  
  }
  Zeresults<-fastReg.func(newY,XX1,AllDNAm)
  return(Zeresults)
}


####################################################################################
#Analyze each self reported ethnicicity
####################################################################################

results.ESS<-myfunc(thepheno$ESS,sqrt)

results.ESS.race<-vector(mode="list",length=3)
for(i in 1:3){
  results.ESS.race[[i]]<-myfunc(thepheno$ESS,sqrt,zerace[[i]])
}



####################################################################################
#Perform the interaction test for the significant probes
####################################################################################

TheQvalues<-qvalue(results.ESS.race[[2]][,5])$qv
ProbesOfInterest<-which(TheQvalues<0.1)

NameOfprobes<-rownames(results.ESS.race[[2]])[ProbesOfInterest]


JustTopProbes<-t(normData.2[as.character(NameOfprobes),])


newcovar<-func.grab.cov(thepheno)
therace  <-as.factor(thepheno[,"race"])
theSite<-as.character(thepheno$Methylation_site)

ESS<-sqrt(thepheno$ESS)

RunTest<-matrix(nrow=ncol(JustTopProbes),ncol=3)

PCs.want<-as.matrix(PC.all)[,1:5]


for(i in 1:ncol(JustTopProbes)){
  Mod1<-lm(ESS~JustTopProbes[,i]+therace+newcovar+theSite+PCs.want)
  E1<-coef(summary(Mod1))[2,c(1,4)]
  Mod2<-lm(ESS~JustTopProbes[,i]*therace+newcovar+theSite+PCs.want)
  E2<-anova(Mod1,Mod2)
  RunTest[i,]<-c(E1,E2[2,"Pr(>F)"])
  
  
}
rownames(RunTest)<-colnames(JustTopProbes)

####################################################################################
#Save out all the results
####################################################################################


save(results.ESS.race,results.ESS,RunTest,TheQvalues,
     file=file.path(results_folder, "data_fin_sub_ESS_betavals_combat_all_includeSNP_PCs_race.Rdata"))



}
q("no")