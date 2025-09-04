

getwd()
install.packages("goseq_1.54.0.tgz", repos = NULL, type = "source")

library(goseq)

# Check the supported genomes
supported_genomes <- supportedGenomes()
supported_genomes$species


# Identify rows containing "Mouse"
mouse_genomes_rows <- supported_genomes$species == "Mouse"

# Print the entire rows for mouse genomes
print(supported_genomes[mouse_genomes_rows, ])





library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
edata = log2(edata+1)


# perform a differential expression analysis using limma
mod = model.matrix(~ pdata_bot$strain)
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma)
limma_pvals = topTable(ebayes_limma,number=dim(edata)[1], adjust.method ="BH", p.value=0.05, sort.by='none')

# first DE gene
limma_pvals[1,]

# number of genes are differentially expressed at the 5% FPR level 
dim(limma_pvals)






library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)



# limma fit with p-value less than 0.05
limma_table = topTable(ebayes_limma,number=dim(edata)[1], adjust.method ="BH", sort.by='none')
genes = as.integer(limma_table$adj.P.Val < 0.05)
names(genes) = rownames(edata)
not_na = !is.na(genes)
genes = genes[not_na]

# use nullp and goseq to perform a gene ontology analysis
pwf = nullp(genes, "mm9", "ensGene")
GO.wall = goseq(pwf, "mm9", "ensGene")
GO.top10 = GO.wall[1:10,1]

# top category
GO.top10[1]

GO.wall$term[1]


library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]



# perform a differential expression analysis using limma, adjusting for lane as a factor
mod_adj = model.matrix(~ pdata_bot$strain + as.factor(pdata_bot$lane.number))
fit_limma_adj = lmFit(edata,mod_adj)
ebayes_limma_adj = eBayes(fit_limma_adj)

# find genes significant at 5% FPR rate
limma_table = topTable(ebayes_limma_adj, number=dim(edata)[1], adjust.method ="BH", sort.by='none')
genes = as.integer(limma_table$adj.P.Val < 0.05)
names(genes) = rownames(edata)
not_na = !is.na(genes)
genes = genes[not_na]

pwf = nullp(genes, "mm9", "ensGene")
GO.wall = goseq(pwf, "mm9", "ensGene")
GO.top10_adj = GO.wall[1:10,1]

# top 10 overrepresented categories are the same
intersect(GO.top10, GO.top10_adj)

## --------------
## [1] "GO:0007129" "GO:0070192" "GO:0045143" "GO:0007127"
