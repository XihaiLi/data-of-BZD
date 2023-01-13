######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("D:\\biowolf\\geoPharm\\18.KEGG")              #���ù���Ŀ¼
rt=read.table("id.txt",sep="\t",header=T,check.names=F)            #��ȡid.txt�ļ�
rt=rt[is.na(rt[,"entrezID"])==F,]                                  #ȥ������idΪNA�Ļ���
gene=rt$entrezID

#kegg��������
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #��������
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)                          #���渻�����

#��״ͼ
pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()

#����ͼ
pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 20,orderBy = "GeneRatio")
dev.off()

######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056