library(methylKit)
library(genomation)
filenames<-list.files(pattern="*.sam")
rufile=as.list(filenames)
ruobj<-processBismarkAln(location=rufile,sample.id=list("Ruskin_Group2","Ruskin_Group1","Ruskin_Group2","Ruskin_Group2","Ruskin_Group2","Ruskin_Group1","Ruskin_Group1","Ruskin_Group1","Ruskin_Group2","Ruskin_Group1","Ruskin_Group1","Ruskin_Group2"),assembly="rn6",
                       save.folder=NULL,save.context=NULL,read.context="CpG",
                       nolap=FALSE,phred64=FALSE,
                       treatment=c(1,0,1,1,1,0,0,0,1,0,0,1))
rumeth=unite(ruobj, destrand=FALSE)
rutiles=tileMethylCounts(ruobj,win.size=1000,step.size=1000)
ruDiff=calculateDiffMeth(rumeth,num.cores=10)
ruDiff25p.hyper=getMethylDiff(ruDiff,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases
ruDiff25p.hypo=getMethylDiff(ruDiff,difference=25,qvalue=0.01,type="hypo")
# get all differentially methylated bases
ruDiff25p=getMethylDiff(ruDiff,difference=25,qvalue=0.01)
gene.obj=readTranscriptFeatures("genes.bed")
rudiffAnn=annotateWithGeneParts(as(ruDiff25p,"GRanges"),gene.obj)
save(list=ls(),file="Ruskin_EAT.rdata",compress=T)
