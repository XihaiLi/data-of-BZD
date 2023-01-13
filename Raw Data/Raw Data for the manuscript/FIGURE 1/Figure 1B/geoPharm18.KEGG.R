######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

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

setwd("D:\\biowolf\\geoPharm\\18.KEGG")              #设置工作目录
rt=read.table("id.txt",sep="\t",header=T,check.names=F)            #读取id.txt文件
rt=rt[is.na(rt[,"entrezID"])==F,]                                  #去除基因id为NA的基因
gene=rt$entrezID

#kegg富集分析
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #富集分析
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
write.table(KEGG,file="KEGG.txt",sep="\t",quote=F,row.names = F)                          #保存富集结果

#柱状图
pdf(file="barplot.pdf",width = 10,height = 7)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()

#气泡图
pdf(file="bubble.pdf",width = 10,height = 7)
dotplot(kk, showCategory = 20,orderBy = "GeneRatio")
dev.off()

######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
