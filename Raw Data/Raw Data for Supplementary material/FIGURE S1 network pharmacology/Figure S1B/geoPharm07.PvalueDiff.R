######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

logFCfilter=0.5                     #logFC������ֵ
pvalueFilter=0.05                   #pֵ��ֵ
conNum=15                           #��������Ʒ��Ŀ
treatNum=15                         #ʵ������Ʒ��Ŀ

library(limma)                    #���ð�
setwd("D:\\biowolf\\geoPharm\\07.diffPvalue")               #���ù���Ŀ¼
rt=read.table("sampleExp.txt",sep="\t",header=T,check.names=F)      #��ȡ�����ļ�

#���һ��������ڶ��У�ȡ��ֵ
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

#ȡlog2����������ļ�����ֵ�ܴ�������Ҫ������ȡlog2
#rt=log2(rt+1)

#�����ݽ���
rt=normalizeBetweenArrays(rt)

#�������
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="all.xls",sep="\t",quote=F)

#���������
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < pvalueFilter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="diff.xls",sep="\t",quote=F,col.names=F)

diffUp <- allDiff[with(allDiff, (logFC>logFCfilter & P.Value < pvalueFilter )), ]
diffUpOut=rbind(id=colnames(diffUp),diffUp)
write.table(diffUpOut,file="up.xls",sep="\t",quote=F,col.names=F)

diffDown <- allDiff[with(allDiff, (logFC < -logFCfilter & P.Value < pvalueFilter )), ]
diffDownOut=rbind(id=colnames(diffDown),diffDown)
write.table(diffDownOut,file="down.xls",sep="\t",quote=F,col.names=F)

#���Ʋ��������ͼ
library(pheatmap)
geneNum=20
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>40){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=rt[hmGene,]
Type=c(rep("C",conNum),rep("T",treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",height=5.5,width=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("green", "black", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         #scale="row",
         fontsize = 10,
         fontsize_row=8,
         fontsize_col=10)
dev.off()

#��ɽͼ
pdf(file="vol.pdf",width=5,height=5)
yMax=max(-log10(allDiff$P.Value))
yMax=ifelse(yMax>100,100,yMax)
xMax=max(abs(allDiff$logFC))
xMax=ifelse(xMax>10,10,xMax)
plot(allDiff$logFC, -log10(allDiff$P.Value), ylab="-log10(P.Value)",xlab="logFC",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1)
diffSub=subset(allDiff, P.Value<pvalueFilter & logFC>logFCfilter)
points(diffSub$logFC, -log10(diffSub$P.Value), pch=20, col="red",cex=1.2)
diffSub=subset(allDiff, P.Value<pvalueFilter & logFC<(-logFCfilter))
points(diffSub$logFC, -log10(diffSub$P.Value), pch=20, col="green",cex=1.2)
abline(v=0,lty=2,lwd=3)
dev.off()

######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056