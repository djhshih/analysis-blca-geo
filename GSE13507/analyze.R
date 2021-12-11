# Purpose: Compare SMARCB1 low vs. SMARC1 high bladder cancers
# from GEO expression data

#---- Preamble ---------------------------------------------------------------
# install_bitbucket("djhshih/mmalign") # install mmalign

library(GEOquery)
library(ggplot2)
library(devtools)
library(mmalign)
library(io)
library(filenamer)
library(limma)
library(dplyr)

accession <- "GSE13507";

target <- "SMARCB1";
probe <- "ILMN_1758823";
cutpoint <- 9.7;

out.fname <- filename(tolower(accession), date=NA, path="out");


#---- Input ------------------------------------------------------------------

gds <- getGEO(accession);
eset <- gds[[1]];

head(pData(phenoData(eset))) # clinical information
head(exprs(eset)) # expression matrix
head(pData(featureData(eset))) 
pheno <- pData(phenoData(eset));

#---- Select only primary cancer samples (165 samples)
pheno_primary_cancer <- pheno[which(pheno$source_name_ch1 == 'Primary bladder cancer'),]
eset_primary_cancer <- eset[,pheno_primary_cancer$geo_accession]

#---- Select only normal samples (9 samples)
pheno_normal <- pheno[which(pheno$source_name_ch1 == 'Normal bladder mucosae'),]
eset_normal <- eset[,pheno_normal$geo_accession]

#---- Select only Recurrent samples (23 samples)
pheno_recurrent <- pheno[which(pheno$source_name_ch1 == 'Recurrent non-muscle invasive tumor'),]
eset_recurrent <- eset[,pheno_recurrent$geo_accession]

#---- PCA analysis to detect potential outliers

pc <- mmalign::pca(exprs(eset_primary_cancer));

qdraw(
  mmalign::pca_plot(pc),
  file = insert(out.fname, tag="pca12", ext="pdf")
)

qdraw(
  mmalign::pca_plot(pc, dims=3:4),
  file = insert(out.fname, tag="pca34", ext="pdf")
)

idx <- ((pc$Z["PC1", ] < -200 & pc$Z["PC2", ] > 150) & (pc$Z["PC3", ] < -100 & pc$Z["PC4", ] > 100));
outliers <- names(which(idx));

# GSM340606 is an outlier which is consistent in PC1&PC2 and PC3&PC4.
#---- Remove GSM340606
pheno_primary_cancer_rm <- pheno_primary_cancer[!(row.names(pheno_primary_cancer) %in% c('GSM340606')),]
eset_primary_cancer_rm <- eset_primary_cancer[,row.names(pheno_primary_cancer)
                                              [-which(row.names(pheno_primary_cancer) == 'GSM340606')]]

#---- Determine which probe of the target gene to use

fidx <- grep(target, pData(featureData(eset))[["Symbol"]]); # get the probe id of SMARCB1 (ILMN_1758823)
pData(featureData(eset))[fidx, ]

target.probes.expr.primary <- t(exprs(eset_primary_cancer_rm)[fidx, ]); # remove the outlier for histogram
target.probes.expr.normal <- t(exprs(eset_normal)[fidx, ]);
target.probes.expr.recurrent <- t(exprs(eset_recurrent)[fidx, ]);

#---- Determine cutpoint
#---- for splitting samples into low and high expression groups

# histogram for primary cancer samples
qdraw(
  {
    hist(target.probes.expr.primary, breaks=40, main=paste0(target,'_primary'), xlab="log expression")
    #abline(v=cutpoint, col="red")
  },
  file = insert(out.fname, tag=c(tolower(target), "probe", "hist_primary"), ext="pdf")
)

# histogram for normal samples
qdraw(
  {
    hist(target.probes.expr.normal, breaks=20, main=paste0(target,'_normal'), xlab="log expression")
    #abline(v=cutpoint, col="red")
  },
  file = insert(out.fname, tag=c(tolower(target), "probe", "hist_normal"), ext="pdf")
)

# histogram for recurrent samples
qdraw(
  {
    hist(target.probes.expr.recurrent, breaks=20, main=paste0(target,'_recurrent'), xlab="log expression")
    #abline(v=cutpoint, col="red")
  },
  file = insert(out.fname, tag=c(tolower(target), "probe", "hist_recurrent"), ext="pdf")
)

# use the lowest point in normal samples as cut point (9.7)
cutpoint <- min(target.probes.expr.normal);

#---- Assign expression groups

groups <- factor(
  as.integer(target.probe.expr >= cutpoint),
  levels=0:1,
  labels=paste0(target, "-", c("low", "high"))
);
pheno_primary_cancer_rm$group <- groups;

qwrite(
  pheno_primary_cancer_rm,
  insert(out.fname, tag=c("clinical_information", tolower(target)), ext="tsv")
);

#---- Assess overlay of expression group on PCA plot

pc_rm <- mmalign::pca(exprs(eset_primary_cancer_rm));
qdraw(
  mmalign::pca_plot(pc_rm, pheno_primary_cancer_rm, aes(colour=group)) +
    coord_fixed()
  ,
  width = 6,
  file = insert(out.fname, tag=c("pca12", "group"), ext="pdf")
)

qdraw(
  mmalign::pca_plot(pc_rm, pheno_primary_cancer_rm, aes(colour=group), dims=3:4) +
    coord_fixed()
  ,
  width = 6,
  file = insert(out.fname, tag=c("pca34", "group"), ext="pdf")
)


#---- Differential expression analysis

design <- data.frame(
  high = 1,
  low_vs_high = ifelse(pheno_primary_cancer_rm$group == "SMARCB1-low", 1, 0)
);

fit <- lmFit(eset_primary_cancer_rm, design);
fit <- eBayes(fit);
res <- topTable(fit, coef="low_vs_high", n=Inf);

res[1:50, c("Symbol", "logFC", "AveExpr", "t", "adj.P.Val")]

# In the SMARCB1-high group,

qwrite(
  select(res,
         ID, Symbol, Entrez_Gene_ID,
         logFC, AveExpr, t, P.Value, adj.P.Val
  ),
  insert(out.fname, tag=c("limma-ebayes", tolower(target), "low_vs_high"), ext="tsv")
);

#---- Prepare input for Broad GSEA

rnk <- res[, c("Symbol", "t")];
idx <- grep("///", rnk$Symbol); # nothing

rnk.c <- rnk;
rnk.c <- rnk.c[order(abs(rnk.c$t), decreasing=TRUE), ];
rnk.c <- rnk.c[rnk.c$Symbol != "", ];
rnk.c <- rnk.c[!duplicated(rnk.c$Symbol), ];
rnk.c <- rnk.c[!is.na(rnk.c$t), ];
#plot(rnk.c$t, pch=".")

qwrite(rnk.c, insert(out.fname, tag=c("limma-ebayes", tolower(target), "low_vs_high"), ext="rnk"));


#---- Correlate group with clinical features

fisher.test(table(pheno_primary_cancer_rm$group, pheno_primary_cancer_rm$"characteristics_ch1.6"))
fisher.test(table(pheno_primary_cancer_rm$group, pheno_primary_cancer_rm$"characteristics_ch1.7"))
fisher.test(table(pheno_primary_cancer_rm$group, pheno_primary_cancer_rm$"characteristics_ch1.8"))
fisher.test(table(pheno_primary_cancer_rm$group, pheno_primary_cancer_rm$"characteristics_ch1.9"))

