library(affy)
setwd("C:/Users/Ling/Desktop/estrogen")
targetsFile <- "estrogen.txt"
pd <- read.AnnotatedDataFrame(targetsFile,header=TRUE,sep="",row.names=1)

ER <- pData(pd)$estrogen
Time <- factor(pData(pd)$time.h)
design <- model.matrix(~ER+Time)
design

design2 <- model.matrix(~ER*Time)
design2

raw <-ReadAffy(celfile.path = "C:/Users/Ling/Desktop/estrogen", filenames=rownames(pData(pd)),phenoData = pd)
raw

eset <- rma(raw)

library(limma)
fit1 <- lmFit(eset, design) # Fit of expression values on a linear line
fit1 <- eBayes(fit1)  # Bayesian algorithm to calculate P values and fold change
topTable(fit1, coef=2)


fit2 <- lmFit(eset, design2)
fit2 <- eBayes(fit2)
topTable(fit2, coef=2)

# Gene(filter)
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")

# Annotation of probe Ids
library(genefilter)

library(GEOquery)
library(limma)
url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33126/matrix/GSE33126_series_matrix.txt.gz"
filenm <- "GSE33126_series_matrix.txt.gz"
if(!file.exists(filenm)) download.file(url, destfile=filenm)
gse <- getGEO(filename=filenm)

# varFilter, fData, anno
gse.expfilter <- varFilter(gse)
anno <- fData(gse.expfilter)
anno <- anno[,c("Symbol", "Entrez_Gene_ID", "Chromosome", "Cytoband")]
fit2$genes <- anno
topTable(fit2)

# Make a volcano plot on fit2 object

volcanoplot(fit2)

