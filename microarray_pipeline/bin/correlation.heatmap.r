args <- commandArgs(TRUE)

#library("RColorBrewer")

Matrix_f = args[1]

mat = read.table(file = Matrix_f, header=T)
dat_matrix = data.matrix(mat[,seq(2,ncol(mat))])
dat_matrix = cor(dat_matrix, method=c("pearson"))
print(dat_matrix)
#dat_matrix = abs(dat_matrix)
#dat_matrix = cor(dat_matrix)

library("gplots")
pdf(file = paste(Matrix_f,".pdf",sep=""))

# plan 1
#heatmap.2(dat_matrix,trace="none", dendrogram="both", col=rev(redblue(100)), margins = c(8, 8), keysize = 1, symkey=FALSE, density.info="none", Rowv=T,Colv=T, cellnote=round(dat_matrix,2), notecol="black", notecex=0.5)
#heatmap.2(dat_matrix,trace="none", dendrogram="both", col=rev(redblue(100)), margins = c(8, 8), keysize = 1, symkey=FALSE, density.info="none", Rowv=T,Colv=T, cellnote=round(dat_matrix,1), notecex=0.2, notecol="black", cexRow = 0.2, cexCol=0.2)

heatmap.2(dat_matrix,trace="none", dendrogram="both", col=rev(redblue(100)), margins = c(8, 8), keysize = 1, symkey=FALSE, density.info="none", Rowv=T,Colv=T, cexRow = 1, cexCol=1, notecol="black")

write.table(dat_matrix,file=paste(Matrix_f,".tab",sep=""), quote=FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
# plan 2
#my_palette <- colorRampPalette(c("blue","white","red"))(n=599)
#col_breaks = c(seq(0,0.33,length=200),seq(0.331,0.66, length=200),seq(0.661,1,length=200))
#heatmap.2(dat_matrix,trace="none", dendrogram="none", col=my_palette, breaks=col_breaks, margins = c(8, 8), keysize = 1, symkey=FALSE, density.info="none", Rowv=F,Colv=F)

# plan X
# cellnote = dat_matrix,
#heatmap.2(dat_matrix,trace="none",hclustfun = hclust,dendrogram="none", col=my_palette, breaks=col_breaks, margins = c(8, 8), keysize = 1.5, symkey=FALSE, density.info="none", Rowv=F,Colv=F)

dev.off()

#col=brewer.pal(11,"RdBu"), 
