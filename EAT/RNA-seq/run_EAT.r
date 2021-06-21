library(DESeq2)
library(ggplot2)
group<-read.delim("EAT_group_cluster.txt",sep="\t",row.names=1)
df<-read.delim("counts.txt",sep="\t",skip=1,row.names=1)
colnames(df)<-sub('\\_.*.bam','',colnames(df))
colnames(df)<-sub('X','',colnames(df))
dd<-df[,6:ncol(df)]
dds<-DESeqDataSetFromMatrix(dd,DataFrame(condition),~condition)
dds<-DESeq(dds,minReplicatesForReplace=Inf)
rld<-rlogTransformation(dds)
pca<-plotPCA(rld,return=T)
percentVar <- round(100 * attr(pca, "percentVar"))
pca$sample<-group$Group
ggplot(pca,aes(PC1,PC2,color=condition,shape=sample))+geom_point(size=4)+xlab(paste0("PC1: ",percentVar[1],"%"))+ylab(paste("PC2: ",percentVar[2],"%"))
ggsave(file="EAT_new.pdf")
######
resbo<-results(dds,contrast = c("condition","Boison_Group1","Boison_Group2"))
reslu<-results(dds,contrast = c("condition","Lubin_Group2","Lubin_Group1"))
resru<-results(dds,contrast = c("condition","Ruskin_Group2","Ruskin_Group1"))
#####
library(richR)
rtgo<-buildAnnot(species="rat",keytype="ENSEMBL",anntype="GO")
rtko<-buildAnnot(species="rat",keytype="ENSEMBL",anntype="KEGG")
###Note we use transcript gene id so we need to replace the gene id with transcript id in the annotation file

######
boko<-richKEGG(rownames(subset(resbo,padj<0.05)),rtko)
bogo<-richGO(rownames(subset(resbo,padj<0.05)),rtgo)
