# Purpose: Compare SMARCB1 low vs. SMARCB1 high bladder cancers
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

accession <- "GSE32894";

target <- "SMARCB1";
cutpoint <- 0.7;

out.fname <- filename(tolower(accession), date=NA, path="out");


#---- Input ------------------------------------------------------------------

gds <- getGEO(accession);
eset <- gds[[1]];

head(pData(phenoData(eset))) # clinical information
head(exprs(eset)) # expression matrix
head(pData(featureData(eset))) 

pheno <- pData(phenoData(eset)); # clinical information

#---- PCA analysis to detect potential outliers (PCA is sensitive to outliers)

pc <- mmalign::pca(exprs(eset));

qdraw(
  mmalign::pca_plot(pc),
  file = insert(out.fname, tag="pca12", ext="pdf")
)

qdraw(
  mmalign::pca_plot(pc, dims=3:4),
  file = insert(out.fname, tag="pca34", ext="pdf")
)

# No obvious outliers.

#---- Determine which probe of the target gene to use

fidx <- grep(target, pData(featureData(eset))[["Symbol"]]); # get the probe id of SMARCB1 (two probes)
pData(featureData(eset))[fidx, ]

target.probes.expr <- t(exprs(eset[fidx, ]));
cor(target.probes.expr) # the correlation level is high

qdraw(
  {
    plot(target.probes.expr)
  },
  file = insert(out.fname, tag=c(tolower(target), "probes", "cor"), ext="pdf") 
);

#---- Take an average of these two probes (ILMN_1758823,ILMN_2403458)
probe_average <- rowMeans(target.probes.expr[,c('ILMN_1758823', 'ILMN_2403458')])

#---- Determine cutpoint
#     for splitting samples into low and high expression groups

target.probe.expr <- probe_average;

qdraw(
  {
    hist(target.probe.expr, breaks=80, main=target, xlab="log expression")
    abline(v=cutpoint, col="red")
  },
  file = insert(out.fname, tag=c(tolower(target), "probe", "hist"), ext="pdf")
)

#---- Assign expression groups

groups <- factor(
  as.integer(target.probe.expr >= cutpoint),
  levels=0:1,
  labels=paste0(target, "-", c("low", "high"))
);
pheno$group <- groups;

qwrite(
  pheno,
  insert(out.fname, tag=c("clinical_information", tolower(target)), ext="tsv")
);

#---- Assess overlay of expression group on PCA plot

qdraw(
  mmalign::pca_plot(pc, pheno, aes(colour=group)) +
    coord_fixed()
  ,
  width = 6,
  file = insert(out.fname, tag=c("pca12", "group"), ext="pdf")
)

qdraw(
  mmalign::pca_plot(pc, pheno, aes(colour=group), dims=3:4) +
    coord_fixed()
  ,
  width = 6,
  file = insert(out.fname, tag=c("pca34", "group"), ext="pdf")
)


#---- Differential expression analysis

design <- data.frame(
  high = 1,
  low_vs_high = ifelse(pheno$group == "SMARCB1-low", 1, 0)
);

fit <- lmFit(eset, design);
fit <- eBayes(fit);
res <- topTable(fit, coef="low_vs_high", n=Inf);

res[1:50, c("Symbol", "logFC", "AveExpr", "t", "adj.P.Val")]

qwrite(
  select(res,
         ID, Symbol, Entrez_Gene_ID,
         logFC, AveExpr, t, P.Value, adj.P.Val
  ),
  insert(out.fname, tag=c("limma-ebayes", tolower(target), "low_vs_high"), ext="tsv")
);

#---- Prepare input for Broad GSEA

rnk <- res[, c("Symbol", "t")];

idx <- grep("///", rnk$Gene.Symbol); # nothing

rnk.c <- rnk
rnk.c <- rnk.c[order(abs(rnk.c$t), decreasing=TRUE), ];
rnk.c <- rnk.c[rnk.c$Symbol != "", ];
rnk.c <- rnk.c[!duplicated(rnk.c$Symbol), ];
rnk.c <- rnk.c[!is.na(rnk.c$t), ];
#plot(rnk.c$t, pch=".")

qwrite(rnk.c, insert(out.fname, tag=c("limma-ebayes", tolower(target), "low_vs_high"), ext="rnk"));

#---- Correlate group with clinical features

fisher.test(table(pheno$group, pheno$"characteristics_ch1.5"))
fisher.test(table(pheno$group, pheno$"characteristics_ch1.6"))

