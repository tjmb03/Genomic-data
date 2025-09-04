
knitr::opts_chunk$set(cache=TRUE)
x = rnorm(10)
plot(x,pch=19,col="dodgerblue")

y = rbinom(20, size=1, prob=0.5)
table(y)









con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

BiocManager::install("plotrix")
library(plotrix)
pie3D(pdata_bm$num.tech.reps,labels=pdata_bm$tissue.type)





con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)

row_sums = rowSums(edata)
index = which(rank(-row_sums) < 500 )
heatmap(edata[index,],Rowv=NA,Colv=NA)


row_sums = rowSums(edata)
edata = edata[order(-row_sums),]
index = 1:500
heatmap(edata[index,],Rowv=NA,Colv=NA)





con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata = pData(bm)
edata = exprs(bm)

# make an MA-plot
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)


library(DESeq2)
rld <- rlog(exprs(bm))

y_rld = rld[,1] - rld[,2]
x_rld = rld[,1] - rld[,2]
plot(x_rld, y_rld, col = "blue", type = "p")


library(DESeq2)
library(ggplot2)

con <- url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file = con)
close(con)

bm <- bodymap.eset
pdata <- pData(bm)
edata <- exprs(bm)

# Log2 transform
log2_edata <- log2(edata + 1)

# rlog transform
dds <- DESeqDataSetFromMatrix(countData = edata, colData = pdata, design = ~ 1)
dds <- DESeq(dds)
rld_edata <- assay(rlog(dds))

# MA-plot with log2 transform
log2_df <- data.frame(A = rowMeans(log2_edata[, c(1, 2)]), M = log2_edata[, 1] - log2_edata[, 2])
ggplot(log2_df, aes(x = A, y = M)) + geom_point() + ggtitle("MA-plot with log2 transform")

# MA-plot with rlog transform
rlog_df <- data.frame(A = rowMeans(rld_edata[, c(1, 2)]), M = rld_edata[, 1] - rld_edata[, 2])
ggplot(rlog_df, aes(x = A, y = M)) + geom_point() + ggtitle("MA-plot with rlog transform")

