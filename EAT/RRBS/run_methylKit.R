###This is the script for Site#2 and Site#3 test
### other two script used same command with different name
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
###read annotation files
gene.obj=readTranscriptFeatures("genes.bed")
###annotate the CpG sites
rudiffAnn=annotateWithGeneParts(as(ruDiff25p,"GRanges"),gene.obj)
save(list=ls(),file="rl_rdata",compress=T)
