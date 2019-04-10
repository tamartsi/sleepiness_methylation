require(GWASTools)
mesa.probes <- getobj("Follow_up_results/20190116_MESA_top_probes.RData")
meta.probes <- getobj("Follow_up_results/20190116_META_AA_top_probes.RData")

all.probes <- c(mesa.probes, meta.probes)


load('Annotation/data_HumanMethylation450_15017582.Rdata')




SubAnnot<-subset(Annotinfo,Name %in% all.probes )[,c("Name","CHR","MAPINFO","UCSC_RefGene_Name")]
#Unique genes
FuncUse<-function(x){
  A1<-unique(strsplit(as.character(x),";")[[1]])
  return(paste(A1,collapse=","))
}

NewGene<-unlist(lapply(SubAnnot$UCSC_RefGene_Name,FuncUse))
SubAnnot<-data.frame(SubAnnot,NewGene)
SubAnnot<-SubAnnot[match(all.probes,SubAnnot$Name),]
all.equal(as.character(SubAnnot$Name),all.probes) #TRUE

top.probes.annot <- data.frame(probeID = as.character(SubAnnot$Name), Chr = SubAnnot$CHR, 
								Position = SubAnnot$MAPINFO, Gene = as.character(SubAnnot$NewGene), stringsAsFactors = FALSE)

save(top.probes.annot, file = "Follow_up_results/20190116_annot_top_probes.RData")