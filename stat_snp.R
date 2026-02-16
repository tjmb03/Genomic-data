library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# recode 0 values to NA
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA

# fit a linear model
lm3 = lm(status ~ snp3)
tidy(lm3)

# fit a logistic regression model
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3)




par(mfrow=c(1,2))

plot(status ~ snp3,pch=19)
abline(lm3,col="darkgrey",lwd=5)
plot(glm3$residuals)





# fit a logistic regression model
snp10 = as.numeric(snpdata[,10])
snp10[snp10==0] = NA
glm10 = glm(status ~ snp10, family="binomial")
tidy(glm10)

snp10_dom = (snp10 == 2)
glm10_dom = glm(status ~ snp10_dom, family="binomial")
tidy(glm10_dom)


# fit an additive logistic regression model to each SNP
results = rep(NA, dim(snpdata)[2])
for (i in 1:ncol(snpdata)){
  snpdata_i = as.numeric(snpdata[,i])
  snpdata_i[snpdata_i == 0] = NA
  glm_i = glm(status ~ snpdata_i, family = "binomial")
  results[i] = tidy(glm_i)$statistic[2]
}

# average effect size
mean(results)

## ---------------
## [1] 0.007155377


# minimum effect size
min(results)

## ---------------
## [1] -4.251469

# maximum effect size
max(results)



# square the coefficients
results_coeff_squre =  results^2

# correlation with the results from using snp.rhs.tests and chi.squared
glm_all = snp.rhs.tests(status ~ 1, snp.data = sub.10)
cor(results_coeff_squre, chi.squared(glm_all))

## ---------------
## [1] 0.9992946



con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)


library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)

# Get the current working directory
current_directory <- getwd()

# Print the current working directory
print(current_directory)



install.packages("genefilter_1.84.0.tgz", repos = NULL, type = "source")
library(genefilter)
install.packages("sva_3.50.0.tgz", repos = NULL, type = "source")
library(sva)
install.packages("edge_2.34.0.tgz", repos = NULL, type = "source")
library(edge)


edata = log2(as.matrix(edata) + 1)

# perform rowttests
tstats_obj = rowttests(edata, as.factor(pdata$population))
tidy(tstats_obj)
# perform rowFtests
fstats_obj = rowFtests(edata, as.factor(pdata$population))
tidy(fstats_obj)
par(mfrow=c(1,2))
hist(tstats_obj$statistic, col=2)
hist(fstats_obj$statistic, col=2)




con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

library(DESeq2)
library(limma)
library(edge)
library(genefilter)

# using DESeq2 test the differences between the studies
de = DESeqDataSetFromMatrix(edata, pdata, ~study)
glm_de = DESeq(de)
result_de = results(glm_de)

result_de

# using limma test the differences
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ as.factor(pdata$study))
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma) 
top = topTable(ebayes_limma,number=dim(edata)[1], sort.by="none")

# correlation in the statistics between two analyses
cor(result_de$stat, top$t)

## -------------
## [1] 0.9278568

# make an MA-plot
y = cbind(result_de$stat, top$t)
limma::plotMA(y)





# DESeq analysis
fp_bh = p.adjust(result_de$pvalue, method="BH")
sum(fp_bh < 0.05)

print(fp_bh[1:10])
## --------
## [1] 1995

# limma analysis
fp_bh = p.adjust(top$P.Value, method="BH")
sum(fp_bh < 0.05)

## --------
## [1] 2807

