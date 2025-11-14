read_count <- read.csv("C/......RNA_data")
metadata <- read.csv("C/...metadata")

dds <- DESeqDataSetFromMatrix(
  countData = read_count,
  colData = metadata,
  design = ~Genotype
)

##Pre-filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

##Perform DE analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("Genotype", "Ndufs4 KO", "WT"))
res <- lfcShrink(dds, coef = "Genotype_WT_vs_Ndufs4.KO", res = res, type = "ashr")

head(res)
write.csv(as.data.frame(res), "DESeq2_results.csv")


-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Differentially Abundance Analysis

library(tidyverse)
library(limma)
library(impute)


###Proteomics
prot_data <- read.csv("Protein_data.csv", row.names = 1, sep = ";", na.strings = c("NaN", "", "NA"))
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



