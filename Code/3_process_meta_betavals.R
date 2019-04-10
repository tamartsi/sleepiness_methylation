require(qvalue)

load("Annotation/data_HumanMethylation450_15017582.Rdata")

load('MESA_CHS_meta_results/20190110_meta_AA_betavals_combat_all.RData')



meta.AA.sig <- meta.AA[which(meta.AA$qval <= 0.05),]


SubAnnot<-subset(Annotinfo,Name %in% rownames(meta.AA.sig))[,c("Name","CHR","MAPINFO","UCSC_RefGene_Name")]
#Unique genes
FuncUse<-function(x){
  A1<-unique(strsplit(as.character(x),";")[[1]])
  return(paste(A1,collapse=","))
}

NewGene<-unlist(lapply(SubAnnot$UCSC_RefGene_Name,FuncUse))
SubAnnot<-data.frame(SubAnnot,NewGene)
SubAnnot<-SubAnnot[match(rownames(meta.AA.sig),SubAnnot$Name),]
all.equal(as.character(SubAnnot$Name),rownames(meta.AA.sig )) #TRUE

meta.AA.sig$Chr <- SubAnnot$CHR
meta.AA.sig$Position <- SubAnnot$MAPINFO
meta.AA.sig$Gene <- SubAnnot$NewGene

save(meta.AA.sig,file = "Follow_up_results/20190110_meta_AA_betavals_combat_all_qval_05.RData")


probes.top.meta <- meta.AA.sig$probeID
save(probes.top.meta, file = "Follow_up_results/20190116_META_AA_top_probes.RData")



#### now target specifically the ones detected by MESA analysis:
meta.AA$qval.mesa <- qvalue(meta.AA$pval.mesa)$qv

mesa.AA.sig <- meta.AA[which(meta.AA$qval.mesa <= 0.05),]


SubAnnot<-subset(Annotinfo,Name %in% rownames(mesa.AA.sig))[,c("Name","CHR","MAPINFO","UCSC_RefGene_Name")]
#Unique genes
FuncUse<-function(x){
  A1<-unique(strsplit(as.character(x),";")[[1]])
  return(paste(A1,collapse=","))
}

NewGene<-unlist(lapply(SubAnnot$UCSC_RefGene_Name,FuncUse))
SubAnnot<-data.frame(SubAnnot,NewGene)
SubAnnot<-SubAnnot[match(rownames(mesa.AA.sig),SubAnnot$Name),]
all.equal(as.character(SubAnnot$Name),rownames(mesa.AA.sig )) #TRUE

mesa.AA.sig$Chr <- SubAnnot$CHR
mesa.AA.sig$Position <- SubAnnot$MAPINFO
mesa.AA.sig$Gene <- SubAnnot$NewGene

save(mesa.AA.sig,file = "Follow_up_results/20190110_meta_AA_betavals_combat_all_mesa_qval_05.RData")

probes.top.mesa <- rownames(mesa.AA.sig$probeID)
save(probes.top.mesa, file = "Follow_up_results/20190116_MESA_top_probes.RData")
