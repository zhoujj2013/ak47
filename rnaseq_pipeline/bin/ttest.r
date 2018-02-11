g1 = read.table(file="./g1.expr")
g2 = read.table(file="./g2.expr")

#g1 = g1+0.01
#g2 = g2+0.01

g1_mean = apply(g1, 1, mean)
g2_mean = apply(g2, 1, mean)

fc = (g2_mean+1)/(g1_mean+1)

logFC = log2(fc)

pv=c()
for(i in 1:length(rownames(g1))){
	r = t.test(g1[i,], g2[i,], alternative = "two.sided")
	pv = append(pv, r$p.value)
}

adjp=p.adjust(pv, method="hochberg")

dat = cbind(g1_mean, g2_mean, fc, logFC, pv, adjp)

write.table(dat, file="ttest.result", row.names=F, col.names=F, sep="\t")

