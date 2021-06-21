filenames<-as.list(list.files(pattern="*.sam"))
group<-read.delim("EAT_group_cluster.txt",sep="\t",row.names=1)
sample.ids<-as.list(group$Group)
library(methylKit)
x<-group$Group
x[x=="Boison_Group1"]=0
x[x=="Boison_Group2"]=1
x[x=="Lubin_Group1"]=2
x[x=="Lubin_Group2"]=3
x[x=="Ruskin_Group1"]=4
x[x=="Ruskin_Group2"]=5
x<-as.numeric(x)
myobj<-processBismarkAln(location=filenames,sample.id=sample.ids,save.folder=NULL,save.context=NULL,read.context="CpG",nolap=FALSE,phred64=FALSE,treatment=x,assembly="rn6")
mymeth<-unite(myobj, destrand=FALSE)
per<-percMethylation(mymeth)
p<-prcomp(t(per))
pca<-as.data.frame(p$x)
library(ggplot2)
pca$group<-sub('\\..*','',rownames(pca))
ggplot(pca,aes(PC1,PC2,color=group))+geom_point()+theme_light()+xlab('')


