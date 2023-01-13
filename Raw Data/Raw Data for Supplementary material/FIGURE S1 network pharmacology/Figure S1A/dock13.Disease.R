######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("venn")


library(venn)              #���ð�
outFile="Disease.txt"      #����ļ�����
setwd("D:\\biowolf\\dock\\13.Disease")    #���ù���Ŀ¼
files=dir()                        #��ȡĿ¼�������ļ�
files=grep("txt$",files,value=T)   #��ȡ.txt��β���ļ�
geneList=list()

#��ȡ����txt�ļ��еĻ�����Ϣ�����浽geneList
for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile==outFile){next}
    rt=read.table(inputFile,header=F)        #��ȡ�����ļ�
    geneNames=as.vector(rt[,1])              #��ȡ��������
    geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
    uniqGene=unique(geneNames)               #����ȡunique
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=uniqGene
    uniqLength=length(uniqGene)
    print(paste(header[1],uniqLength,sep=" "))
}

#����vennͼ
mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#91612D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="venn.pdf",width=5,height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()

#���沢������
unionGenes=Reduce(union,geneList)
write.table(file=outFile,unionGenes,sep="\t",quote=F,col.names=F,row.names=F)


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056