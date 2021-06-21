######list all files 
filenames<-as.list(list.files(pattern="*.sam"))
##### read group information
group<-read.delim("EAT_group_cluster.txt",sep="\t",row.names=1)
sample.ids<-as.list(group$Group)
library(methylKit)
x<-group$Group
x[x=="Boison_Group2"]=0
x[x=="Lubin_Group1"]=1
x[x=="Ruskin_Group1"]=2
x<-as.numeric(x)
myobj<-processBismarkAln(location=filenames,sample.id=sample.ids,save.folder=NULL,save.context=NULL,read.context="CpG",nolap=FALSE,phred64=FALSE,treatment=x,assembly="rn6")
mymeth<-unite(myobj, destrand=FALSE)
per<-percMethylation(mymeth)
####use all CpG sites or use top 10
p<-prcomp(t(per))
explain <- p$sdev^2/sum(p$sdev^2)
pca<-as.data.frame(p$x)
library(ggplot2)
pca$group<-sub('\\..*','',rownames(pca))
ggplot(pca,aes(PC1,PC2,color=group))+geom_point()+theme_light()+xlab(paste0("PC1: ("explain[1]*100,"%)"))+ylab(paste0("PC2: ("explain[1]*100,"%)"))

