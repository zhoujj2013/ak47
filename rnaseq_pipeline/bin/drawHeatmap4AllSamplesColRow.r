args <- commandArgs(TRUE)

if (length(args) < 4) {
	stop(paste("Rscript XXX.R prefix colnum[1,2,3,4] namestr[a,b,c,d] expr.tab",sep=""))
}

# input max genes number 8000
prefix = args[1]
colnum = args[2]
namestr = args[3]
expr = args[4]

colarr = unlist(strsplit(colnum,"[,]"))
names = unlist(strsplit(namestr,"[,]"))

print("## deal with dataset")
d1 <- read.table(file=expr, header=F)

dat = data.frame()
for (i in 1:length(colarr)) {
	print(i)
	if (i == 1){
		dat = d1[,strtoi(colarr[i])]
	} else
		dat = cbind(dat, d1[,strtoi(colarr[i])])
}
m1 = as.matrix(dat)

#m1<- as.matrix( d1[,2:ncol(d1)])
rownames(m1) = d1$V1
colnames(m1) = names

print("# Sample names")

print(names)

print("# Data summary")
summary(m1[,1])
#m1<- log10(m1)
#summary(m1[,1])
min(m1)
max(m1)

dat_matrix = m1

library("gplots")
hmcols <- rev(redblue(2750));
pdf(paste(prefix,".all.sample.heatmap.pdf",sep=""))
heatmap.2(dat_matrix,trace="none",hclustfun = hclust,dendrogram="both", margins = c(8, 8), keysize = 1.5, symkey=FALSE, density.info="none", Rowv=T,Colv=T, scale="row", col=hmcols)
dev.off()
