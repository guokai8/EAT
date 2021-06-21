###This is the script for Site#2 and Site#3 test
### other two script used same command with different group name
##load packages
library(methylKit)
library(genomation)
##get filenames
filenames<-list.files(pattern="*.sam")
rufile=as.list(filenames)
##read sam files
ruobj<-processBismarkAln(location=rufile,sample.id=list("Lubin_Group1","Ruskin_Group1","Lubin_Group1","Lubin_Group1", "Ruskin_Group1",
                                                        "Lubin_Group1","Ruskin_Group1", "Ruskin_Group1","Lubin_Group1","Ruskin_Group1",
                                                        "Lubin_Group1","Ruskin_Group1")
                       save.folder=NULL,save.context=NULL,read.context="CpG",
                        mincov = 10,minqual = 20,
                       nolap=FALSE,phred64=FALSE,
                       treatment=c(0,1,0,0,1,0,1,1,0,1,0,1))
###merged files
rumeth=unite(ruobj, destrand=FALSE)
clusterSamples(rumeth)
###windows for DMR
rutiles=tileMethylCounts(ruobj,win.size=1000,step.size=1000)
###calculate the differential CpGs
ruDiff=calculateDiffMeth(rumeth,num.cores=30)
# get all differentially methylated bases
ruDiff25p=getMethylDiff(ruDiff,difference=25,qvalue=0.01)
# we just focus on the no sex chromosome
ruDiff25p=ruDiff25p[ruDiff25p$chr%in%c(1:20),]
###read annotation files
gene.obj=readTranscriptFeatures("genes.bed")
###annotate the CpG sites
rudiffAnn=annotateWithGeneParts(as(ruDiff25p,"GRanges"),gene.obj)
rugene<-getAssociationWithTSS(rudiffAnn)
library(richR)
rtgo<-buildAnnot(species="rat",keytype="ENSEMBL",anntype="GO")
rtko<-buildAnnot(species="rat",keytype="ENSEMBL",anntype="KEGG")
#####
rtgo<-as.data.frame(rtgo)
rtko<-as.data.frame(rtko)
###Note we use transcript gene id so we need to replace the gene id with transcript id in the annotation file
geneid<-read.table('EATdata/rat.txt',sep="\t",skip=1)
colnames(geneid)<-c("TranID","GeneID")
library(tidyverse)
rtgo<-left_join(rtgo,geneid,by=c('GeneID'='GeneID'))%>%select(TranID,GOALL,ONTOLOGYAL,Annot)
rtko<-left_join(rtko,geneid,by=c('GeneID'='GeneID'))%>%select(TranID,PATH,Annot)
######
ruko<-richKEGG(unique(rugene$feature.name),rtko)
rugo<-richGO(unique(rugene$feature.name)),rtgo)
write.table(rugo,file="Site2vsSite3_GO.csv")
write.table(ruko,file="Site2vsSite3_KEGG.csv")
save(list=ls(),file="rl_rdata",compress=T)
