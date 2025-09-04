# see vignette for suggestions on generating
# count tables from RNA-Seq data
cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
cond <- factor(rep(1:2, each=5))
# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
# standard analysis
dds <- DESeq(dds)
res <- results(dds)
# moderated log2 fold changes
resultsNames(dds)
resLFC <- lfcShrink(dds, coef=2, type="apeglm")
# an alternate analysis: likelihood ratio test
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)