######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

logFCfilter=0.5                     #logFC过滤阈值
pvalueFilter=0.05                   #p值阈值
conNum=15                           #正常组样品数目
treatNum=15                         #实验组样品数目

library(limma)                    #引用包
setwd("D:\\biowolf\\geoPharm\\07.diffPvalue")               #设置工作目录
rt=read.table("sampleExp.txt",sep="\t",header=T,check.names=F)      #读取输入文件

#如果一个基因存在多行，取均值
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

#取log2，如果输入文件的数值很大，我们需要对数据取log2
#rt=log2(rt+1)

#对数据矫正
rt=normalizeBetweenArrays(rt)

#差异分析
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="all.xls",sep="\t",quote=F)

#输出差异结果
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < pvalueFilter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="diff.xls",sep="\t",quote=F,col.names=F)

diffUp <- allDiff[with(allDiff, (logFC>logFCfilter & P.Value < pvalueFilter )), ]
diffUpOut=rbind(id=colnames(diffUp),diffUp)
write.table(diffUpOut,file="up.xls",sep="\t",quote=F,col.names=F)

diffDown <- allDiff[with(allDiff, (logFC < -logFCfilter & P.Value < pvalueFilter )), ]
diffDownOut=rbind(id=colnames(diffDown),diffDown)
write.table(diffDownOut,file="down.xls",sep="\t",quote=F,col.names=F)

#绘制差异基因热图
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

#火山图
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
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
