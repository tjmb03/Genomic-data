if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALL")


library(ALL)
data(ALL)
mean(exprs(ALL[,5]))

library(biomaRt)

BiocManager::install("hgu95av2.db")
library("hgu95av2.db")
packageDescription("hgu95av2.db")$Version
ls("package:hgu95av2.db")

keytypes(hgu95av2.db)
keys(hgu95av2.db, keytype = "SYMBOL")
columns(hgu95av2.db)
select(hgu95av2.db, keys = c("CMA1", "CMd18"), keytype = "SYMBOL", columns = c("ENTREZID", "ENSEMBL"))
keys(hgu95av2.db, keytype = "PROBEID")






# Load the hgu95av2.db package
library(hgu95av2.db)

# Get all keys for the specified keytype ("affy_hg_u95av2")
all_keys <- keys(hgu95av2.db, keytype = "PROBEID")
all_keys1 <- keys(hgu95av2.db, keytype = "SYMBOL")
# Specify the columns you want to retrieve (replace with actual column names)
desired_columns <- c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME")
desired_columns <- "ENSEMBL"
# Use the select function to retrieve attributes for all keys
result <- select(hgu95av2.db, keys = all_keys, keytype = "PROBEID", columns = desired_columns)

# Print the result
head(result)



library(AnnotationDbi)
# list Marts
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
mart <- useMart(host = 'https://useast.ensembl.org/index.html', biomart = 'ENSEMBL_MART_ENSEMBL')
mart
head(listDatasets(mart))
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
ensembl
# list Attributes
attributes <- listAttributes(ensembl)
head(attributes$name)
filters <- listFilters(ensembl)
head(filters)
# annotate each feature
feature_name <- featureNames(ALL)
head(feature_name)
annotation_ALL <- getBM(attributes = c("ensembl_gene_id","affy_hg_u95av2"), 
                        filters = "affy_hg_u95av2", 
                        values = feature_name, mart = ensembl)


sum(table(annotation_ALL[,2])>1)






# Load required libraries
library(biomaRt)
head(listMarts())


# list Attributes
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

# annotate autosomes
chrom <- c(1:22)
annotation_ALL_chr <- getBM(attributes=c("ensembl_gene_id", "affy_hg_u95av2", "chromosome_name"), 
                            filters=c("affy_hg_u95av2","chromosome_name"), 
                            values=list(feature_name, chrom), 
                            mart=ensembl)

sum(table(table(annotation_ALL_chr[,2])))








# Install the minfiData package
install.packages("minfiData")
library(minfiData)
install.packages("minfi")
library(minfi)


getOption("repos")
pkg <- "minfi"  # Replace "minfi" with the package name you are interested in
av <- available.packages(filters = list())
av[av[, "Package"] == pkg, ]

pkg_info <- available.packages()["minfi", ]
print(pkg_info)

if ("minfi" %in% rownames(installed.packages())) {
  pkg_info <- as.data.frame(installed.packages()["minfi", ])
  print(pkg_info)
} else {
  print("minfi is not installed.")
}





mean(getMeth(MsetEx)[,2])




library(GEOquery)
eList <- getGEO("GSE788")
eList
eData <- eList[[1]]
eData

mean(exprs(eData)[,2])


BiocManager::install("airway")
library(airway)
library(GenomicRanges)
data(airway)
mean(airway$avgLength)


sum(assay(airway)[,3]>=1)




library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

TxDb.Hsapiens.UCSC.hg19.knownGene

# exon data of txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb
txdb_exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_exons
# transcripts on the autosome
autosome <- paste0("chr", c(1:22))
txdb_exons_autosome <- keepSeqlevels(txdb_exons, autosome, pruning.mode = "coarse")

# rename in NCBI format
txdb_ncbi <- mapSeqlevels(seqlevels(txdb_exons), "NCBI")
txdb_exons_ncbi <- renameSeqlevels(txdb_exons_autosome, txdb_ncbi)

dim(subsetByOverlaps(airway, txdb_exons_ncbi))[1]





sample_SRR1039508 <- airway[, 1]
sample_SRR1039508_autosome <- subsetByOverlaps(sample_SRR1039508, txdb_exons_ncbi)

autosome_reads <- sum(assay(sample_SRR1039508_autosome, "counts"))
total_reads <- sum(assay(sample_SRR1039508, "counts"))

# percentage of the total reads in airway dataset for SRR1039508 which overlaps autosome of txdb
autosome_reads/total_reads



library(AnnotationHub)
ah <- AnnotationHub()


ah_E096 <- query(ah, c("E096", "H3K4me3", "narrowPeak"))
ah_record <- ah_E096[["AH30596"]]

ah_record_autosome <- keepSeqlevels(ah_record, autosome, pruning.mode = "coarse")
ah_record_ncbi <- renameSeqlevels(ah_record_autosome, txdb_ncbi)

ncbi_group <- extractSeqlevelsByGroup(species = "Homo sapiens", style = "NCBI", group = "auto")
sample_ncbi <- keepSeqlevels(range(rowRanges(sample_SRR1039508_autosome)), ncbi_group)

ov <- subsetByOverlaps(promoters(sample_ncbi), ah_record_ncbi)
ov <- subsetByOverlaps(sample_SRR1039508, ov)

median(assay(ov, "counts"))
