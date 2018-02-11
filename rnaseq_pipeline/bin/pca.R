args <- commandArgs(TRUE)

color_str = args[1]
colarr = unlist(strsplit(color_str,"[,]"))

library("ggplot2")
data=read.table("degs.expr",header=T,sep="\t",row.names=1)
data.pca <- prcomp(data)
pca.sum = summary(data.pca)

pdf(file="pcaBarplot.pdf",width=15)
barplot(pca.sum$importance[2,1:8]*100,xlab="PC",ylab="percent",col="skyblue")
dev.off()

pca.sum$importance

pdf(file="PCA2d.pdf")
ggplot(data = data.frame(data.pca$rotation), aes(PC1, PC2)) + geom_point(color = colarr)+
annotate("text",x=data.frame(data.pca$rotation)$PC1-0.01,y=data.frame(data.pca$rotation)$PC2,label=rownames(data.pca$rotation))
#geom_path(data=data.frame(data.pca$rotation), aes(x=PC1, y=PC2,color=c("black", "black","black", "black","black", "black","black", "black","red", "red", "red", "red","red", "red", "red", "red","red", "red", "red")), size=1, linetype=2)
dev.off()
