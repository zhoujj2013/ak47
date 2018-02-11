args = commandArgs(TRUE)

category_f=args[1]
outfig=args[2]

library(ggplot2)
dat = read.table(file=category_f, header=F, sep="\t")
dat$V2 = -log10(dat$V2)

#colnames(dat) = c("Category","-log10(pvalue)")

dat$V1 <- factor(dat$V1, levels = dat$V1[order(dat$V2)])

pdf(outfig, width=10,height=4)
p<-ggplot(data=dat, aes(x=V1, y=V2)) +
  geom_bar(stat="identity") + 
  ylab("-log10(p-value)") +
  xlab("Category") +
  coord_flip()
p
dev.off()

