if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

# Install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("AnnotationHub", "Biostrings", "BSgenome", "GenomicRanges", "rtracklayer", "Rsamtools", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene"))



library(AnnotationHub)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Load the genome information
genome <- BSgenome.Hsapiens.UCSC.hg19


# count total bases on chr22
alphabet_frequency <- alphabetFrequency(Hsapiens$chr22)
total_bases <- sum(alphabet_frequency[c('A','G','C','T')])
# count GC bases on chr22  
GC_bases <- sum(alphabet_frequency[c('G','C')])

#calculate GC ratio
GC_content <- GC_bases/total_bases
GC_content


# retrieve record
ah <- AnnotationHub()
H3K27me3_qh <- query(ah, c("H3K27me3", "E003", "narrowPeak"))
H3K27me3_record <- H3K27me3_qh[["AH29892"]]

# extract chr 22 and sequences
H3K27me3_chr22 <- subset(H3K27me3_record, seqnames == "chr22")
H3K27me3_chr22_views <- Views(Hsapiens, H3K27me3_chr22)

# calculate mean GC content
GC_contents <- letterFrequency(H3K27me3_chr22_views, "GC", as.prob = TRUE)
mean_GC <- mean(GC_contents)

mean_GC

signal_value <- mcols(H3K27me3_chr22_views)$signalValue
cor(signal_value, GC_contents)