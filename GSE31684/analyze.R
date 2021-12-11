# Purpose: Compare SMARCB1 low vs. SMARC1 high bladder cancers
# from GEO expression data

#---- Preamble ---------------------------------------------------------------

library(GEOquery)
library(ggplot2)
library(mmalign)
library(io)
library(filenamer)
library(limma)

accession <- "GSE31684";

target <- "SMARCB1";
probe <- "212167_s_at";
cutpoint <- 8.3;

out.fname <- filename(tolower(accession), date=NA, path="out");


#---- Input ------------------------------------------------------------------

gds <- getGEO(accession);
eset <- gds[[1]];

head(pData(phenoData(eset)))
head(exprs(eset))
head(pData(featureData(eset)))

pheno <- pData(phenoData(eset));

#---- PCA analysis to detect potential outliers

pc <- mmalign::pca(exprs(eset));

qdraw(
	mmalign::pca_plot(pc),
	file = insert(out.fname, tag="pca12", ext="pdf")
)

qdraw(
	mmalign::pca_plot(pc, dims=3:4),
	file = insert(out.fname, tag="pca34", ext="pdf")
)

idx <- (pc$Z["PC3", ] < 0 & pc$Z["PC4", ] < -100) | (pc$Z["PC3", ] < -150);
outliers <- names(which(idx));
# For now, we don't remove these putative outliers.

#---- Determine which probe of the target gene to use

fidx <- grep(target, pData(featureData(eset))[["Gene Symbol"]]);
pData(featureData(eset))[fidx, ]

target.probes.expr <- t(exprs(eset[fidx, ]));
cor(target.probes.expr)

qdraw(
	{
		plot(target.probes.expr)
	},
	file = insert(out.fname, tag=c(tolower(target), "probes", "cor"), ext="pdf")
);

#---- Determine cutpoint
#     for splitting samples into low and high expression groups

target.probe.expr <- target.probes.expr[, probe];

qdraw(
	{
		hist(target.probe.expr, breaks=30, main=target, xlab="log expression")
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
	file = insert(out.fname, tag="pca34", ext="pdf")
)


#---- Differential expression analysis

design <- data.frame(
	high = 1,
	low_vs_high = ifelse(pheno$group == "SMARCB1-low", 1, 0)
);

fit <- lmFit(eset, design);
fit <- eBayes(fit);
res <- topTable(fit, coef="low_vs_high", n=Inf);

res[1:50, c("Gene.Symbol", "logFC", "AveExpr", "t", "adj.P.Val")]

# In the SMARCB1-high group,
# many SWI/SNF complex members are upregulated
# many hnRNP complex members are upregulated

qwrite(
	select(res,
		ID, Gene.Symbol, ENTREZ_GENE_ID,
		logFC, AveExpr, t, P.Value, adj.P.Val
	),
	insert(out.fname, tag=c("limma-ebayes", tolower(target), "low_vs_high"), ext="tsv")
);

#---- Prepare input for Broad GSEA

rnk <- res[, c("Gene.Symbol", "t")];

idx <- grep("///", rnk$Gene.Symbol);

# entries with single genes
rnk.sing <- rnk[-idx, ];

# split entries with multiple genes
rnk.mult <- rnk[idx, ];
ss <- strsplit(rnk.mult$Gene.Symbol, " /// ");
rnk.mult1 <- data.frame(
	Gene.Symbol = unlist(lapply(ss, function(x) x[1])),
	t = rnk.mult$t
);
rnk.mult2 <- data.frame(
	Gene.Symbol = unlist(lapply(ss, function(x) x[2])),
	t = rnk.mult$t
);

rnk.c <- rbind(rnk.sing, rnk.mult1, rnk.mult2);
rnk.c <- rnk.c[order(abs(rnk.c$t), decreasing=TRUE), ];
rnk.c <- rnk.c[rnk.c$Gene.Symbol != "", ];
rnk.c <- rnk.c[!duplicated(rnk.c$Gene.Symbol), ];
rnk.c <- rnk.c[!is.na(rnk.c$t), ];
#plot(rnk.c$t, pch=".")

qwrite(rnk.c, insert(out.fname, tag=c("limma-ebayes", tolower(target), "low_vs_high"), ext="rnk"));


#---- Correlate group with clinical features

fisher.test(table(pheno$group, pheno$"rc_stage:ch1"))
fisher.test(table(pheno$group, pheno$"rc grade:ch1"))
fisher.test(table(pheno$group, pheno$"metastasis:ch1"))
fisher.test(table(pheno$group, pheno$"cluster:ch1"))

