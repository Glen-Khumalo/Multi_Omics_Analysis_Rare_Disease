cr <- read.csv("C:/Users/Glen Khumalo/OneDrive - North-West University/Documents/Postdoctoral_2025/Multi_omics intergration/OB_RNASeq_data.csv", header = TRUE, 
               sep =",", row.names = 1)
colnames(cr) <- sub("^X", "", colnames(cr))
head(cr)
mr <- read.csv("C:/Users/Glen Khumalo/OneDrive - North-West University/Documents/Postdoctoral_2025/Multi_omics intergration/OB_RNASeq_metadata.csv", header = TRUE, stringsAsFactors = FALSE,
               sep= ";", row.names = 1)


dds <- DESeqDataSetFromMatrix(
  countData = cr,
  colData = mr,
  design = ~Genotype
)

##Pre-filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

##Perform DE analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("Genotype", "Ndufs4 KO", "WT"))
res <- lfcShrink(dds, coef = "Genotype_WT_vs_Ndufs4.KO", res = res, type = "ashr")

head(res)
write.csv(as.data.frame(res), "DESeq2_results_OB.csv")


-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Differentially Abundance Analysis

library(tidyverse)
library(limma)
library(impute)


## Set wowrking directory
setwd("C:/Users/Glen Khumalo/OneDrive - North-West University/Documents/Postdoctoral_2025/Multi_omics intergration")

###Proteomics
prot_data <- read.csv("OB_SWATH_data.csv", row.names = 1, sep = ";", na.strings = c("NaN", "", "NA"))
prot_data <- prot_data[, !grepl("^X", colnames(prot_data)), drop = TRUE]


# Remove empty columns and columns starting with X
empty_cols <- sapply(prot_data, function(x) all(is.na(x)))
data <- prot_data[, !empty_cols]

# Extract protein IDs and expression matrix
protein_ids <- data$Uniprot_ID
expression_matrix <- as.matrix(data[, -1])
rownames(expression_matrix) <- protein_ids

#Remove proteins with missing values (more tha  50% in a ny group)
wt_cols <- grep("^WT", colnames(expression_matrix))
ko_cols <- grep("^KO", colnames(expression_matrix))

missing_threshold <- 0.5
keep_proteins <- apply(expression_matrix, 1, function(x) {
  wt_missing <- sum(is.na(x[wt_cols])) / length(wt_cols)
  ko_missing <- sum(is.na(x[ko_cols])) / length(ko_cols)
  wt_missing < missing_threshold & ko_missing < missing_threshold
})

filtered_matrix <- expression_matrix[keep_proteins, ]

# Impute missing values using k-nearest neighbors
if(any(is.na(filtered_matrix))) {
  set.seed(123)
  imputed_matrix <- impute.knn(filtered_matrix)$data
} else {
  imputed_matrix <- filtered_matrix
}

filtered_matrix <- as.matrix(filtered_matrix)
storage.mode(filtered_matrix) <- "numeric"
imputed <- impute.knn(filtered_matrix)

# Log2 transformation (common for proteomics data)
log2_matrix <- log2(filtered_matrix + 1)

# Create experimental design
group <- factor(c(rep("WT", length(wt_cols)), rep("KO", length(ko_cols))))
design <- model.matrix(~ 0 + group)
colnames(design) <- c("KO", "WT")

# LIMMA differential expression analysis
fit <- lmFit(log2_matrix, design)

# Create contrast (WT vs KO)
contrast_matrix <- makeContrasts(WTvsKO = WT - KO, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, number = Inf, adjust.method = "fdr")
results$Protein <- rownames(results)


#Remove empty columns
# Remove empty columns and columns starting with X
empty_cols <- sapply(prot_data, function(x) all(is.na(x)))
prot_data <- prot_data[, !empty_cols]

prot_data <- prot_data[, !grepl("^X", colnames(prot_data))]
# Inspect
dim(prot_data)
head(prot_data[, 1:5])

met_prot <- read.csv("OB_SWATH_metadata.csv", row.names = 1, sep = ";")
# Create design matrix
design <- model.matrix(~0 + met_prot$Genotype, data = met_prot)
colnames(design) <- gsub("Genotype", "", colnames(design))
design


# Fit the linear model
fit <- lmFit(prot_data, design)

# Define contrasts (KO vs WT)
contrast.matrix <- makeContrasts(KOvsWT = KO - WT, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)


dir.create("C:/Rlibs", showWarnings = FALSE)
.libPaths("C:/Rlibs")
install.packages("xfun")
BiocManager::install("impute")