######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("venn")


library(venn)              #引用包
outFile="Disease.txt"      #输出文件名称
setwd("D:\\biowolf\\dock\\13.Disease")    #设置工作目录
files=dir()                        #获取目录下所有文件
files=grep("txt$",files,value=T)   #提取.txt结尾的文件
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile==outFile){next}
    rt=read.table(inputFile,header=F)        #读取输入文件
    geneNames=as.vector(rt[,1])              #提取基因名称
    geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
    uniqGene=unique(geneNames)               #基因取unique
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=uniqGene
    uniqLength=length(uniqGene)
    print(paste(header[1],uniqLength,sep=" "))
}

#绘制venn图
mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#91612D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="venn.pdf",width=5,height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()

#保存并集基因
unionGenes=Reduce(union,geneList)
write.table(file=outFile,unionGenes,sep="\t",quote=F,col.names=F,row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056
