library(Biobase)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# No transformations
svd1 = svd(edata)
ori_pca = svd1$d^2/sum(svd1$d^2)
ori_pca[1]
# log2 transform
edata_log2 = log2(edata + 1)
svd2 = svd(edata_log2)
log2_pca = svd2$d^2/sum(svd2$d^2)
log2_pca[1]
# log2 transform, subtract row means
edata_centered = edata_log2 - rowMeans(edata_log2)
svd3 = svd(edata_centered)
centered_data_pca = svd3$d^2/sum(svd3$d^2)
centered_data_pca[1]




# Assuming svd3$d contains the singular values
variance_explained <- svd3$d^2 / sum(svd3$d^2)

# Plot the scree plot
plot(1:length(variance_explained), variance_explained, type = "b", 
     xlab = "Principal Component", ylab = "Variance Explained",
     main = "Scree Plot")

# Add labels and grid
abline(h = 0.05, col = "red", lty = 2)  # Add a threshold line if needed
grid()


# use svd to calculate the singular vectors
set.seed(333)
svd1 = svd(edata_centered)

edata_kmeans = kmeans(t(edata_centered), centers=2)
cor.test(svd1$v[,1], edata_kmeans$cluster)


# Assuming svd1 is the result of svd(edata_centered)
# and edata_kmeans is the result of kmeans(t(edata_centered), centers = 2)

# Scatter plot
plot(svd1$v[, 1], edata_kmeans$cluster, col = c("red", "blue")[as.numeric(edata_kmeans$cluster)], pch = 16, xlab = "Singular Vector 1", ylab = "Cluster Assignment")
# Check the levels of edata_kmeans$cluster
levels(edata_kmeans$cluster)
# Convert to factor
edata_kmeans$cluster <- as.factor(edata_kmeans$cluster)

# Check the levels
levels(edata_kmeans$cluster)
# Add a legend
legend("topright", legend = levels(edata_kmeans$cluster), col = c("red", "blue"), pch = 16, title = "Cluster")

# Add title and labels
title(main = "Scatter Plot of Singular Vector 1 vs Cluster Assignment")









con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# fit linear model
lm1 = lm(edata[1,] ~ pdata_bm$num.tech.reps)

table(pdata_bm$num.tech.reps)
pdata_bm$num.tech.reps


# plot the data
plot(pdata_bm$num.tech.reps,edata[1,])
abline(lm1$coeff[1], lm1$coeff[2], col=2, lwd=3)
lm1$coeff[1]
lm1$coeff[2]


# fit linear model
lm2 = lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
summary(lm2)

plot(pdata_bm$age,edata[1,])
abline(lm2$coeff[1], lm2$coeff[2], col=2, lwd=3)
plot(pdata_bm$gender,edata[1,])
abline(lm2$coeff[1], lm2$coeff[2], col=2, lwd=3)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)





edata = log2(edata + 1)

# fit a regression model to each sample, using population as the outcome
mod = model.matrix(~ pdata$population)
fit = lm.fit(mod, t(edata))

# dimension of the residual matrix
dim(fit$residuals)
fr = fit$residuals
# dimension of the effects matrix
dim(fit$effects)
fe = fit$effects
# dimension of the coefficients matrix
dim(fit$coefficients)
fc = fit$coefficients
fv = fit$fitted.values

?lm.fit
head(fit$effects)


con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)


BiocManager::install(c("devtools","Biobase","limma","edgeR",))
library(devtools)
library(Biobase)
library(limma)
library(edge)

# subset the expression data to the samples without mimssing values of age
pdata_bm = na.omit(pdata_bm)
edata = edata[,rownames(pdata_bm), drop=FALSE]

# fit many regression models to the expression data where age is the outcome
mod_adj = model.matrix(~ pdata_bm$age)
fit_limma = lmFit(edata,mod_adj)

fit_limma$coefficients[1000,]
flc = fit_limma$coefficients

# make a plot of the 1,000th gene and fitted values
intercept = fit_limma$coefficients[1000,][1]
slope = fit_limma$coefficients[1000,][2]
x = edata[1000,]*slope+intercept

edata[1000,]

plot(x,pdata_bm$age)





con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

pdata_bm$tissue.type

# subset the expression data to the samples without mimssing values of age
pdata_bm = na.omit(pdata_bm)
edata = edata[,rownames(pdata_bm), drop=FALSE]

pdata_bm$tissue.type


con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

head(pdata)

##         sample.id num.tech.reps population      study
## NA06985   NA06985             1        CEU Montgomery
## NA06986   NA06986             1        CEU Montgomery
## NA06994   NA06994             1        CEU Montgomery
## NA07000   NA07000             1        CEU Montgomery
## NA07037   NA07037             1        CEU Montgomery
## NA07051   NA07051             1        CEU Montgomery





con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)



library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)



# preprocessing the data
set.seed(33353)
pheno = na.omit(pdata_bm)
edata = edata[,rownames(pheno), drop=FALSE]
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 1,]

# fit a sva model
mod = model.matrix(~age, data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva1 = sva(edata, mod,mod0, n.sv=2)

# correlation between surrogate for batch and age
cor(sva1$sv, pheno$age)

# correlation between surrogate for batch and race
cor(sva1$sv, as.numeric(pheno$race))

# correlation between surrogate for batch and gender
cor(sva1$sv, as.numeric(pheno$gender))


browseVignettes("sva")










