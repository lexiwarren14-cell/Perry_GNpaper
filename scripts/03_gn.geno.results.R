#libraries and functions
# setwd("")

library(rrBLUP)
library(dplyr)

source("functions.R")

#read in data
geno <- read.table("../BarleyMarkersConsensusWithChecks2024.txt")
phenoDataAll <- read.csv("../phenoDataAll.csv", check.names = FALSE)

#make year-specific subsets
gn23 <- phenoDataAll[phenoDataAll$Year == 2023, ]
gn24 <- phenoDataAll[phenoDataAll$Year == 2024, ]

gn23$GN <- as.numeric(gn23$GN)
gn23 <- gn23[!is.na(gn23$GN), ]
gn24$GN <- as.numeric(gn24$GN)
gn24 <- gn24[!is.na(gn24$GN), ]
gn24 <- gn24[!is.na(gn24$moisture), ]
gn24 <- gn24[!is.na(gn24$fineExtract), ]
gn24 <- gn24[!is.na(gn24$SP), ]
gn24 <- gn24[!is.na(gn24$FAN), ]
gn24 <- gn24[!is.na(gn24$BG), ]
gn24 <- gn24[!is.na(gn24$AA), ]
gn24 <- gn24[!is.na(gn24$DP), ]

#recombine the filtered year-specific subsets in the same 2023-then-2024 order
pheno23.24 <- rbind(gn23, gn24)

#kinship matrix
K <- make.K(geno)
#checking lines are all there
all(colnames(K) %in% pheno23.24$Line)
all(pheno23.24$Line %in% colnames(K))
#missing some genotypes so getting rid of them
no.geno <- unique(pheno23.24$Line[!pheno23.24$Line %in% colnames(K)])
pheno23.24.2 <- pheno23.24[!pheno23.24$Line %in% no.geno, ]
#making sure K is in the correct order
pheno23.24.2$Line <- factor(pheno23.24.2$Line, levels = unique(pheno23.24.2$Line))
lines <- levels(pheno23.24.2$Line)
K.gwas2 <- K[lines, lines]
all(lines == colnames(K.gwas2))
colnames(K.gwas2) <- cleanBarleyNames(T3ize(colnames(K.gwas2)))
rownames(K.gwas2) <- cleanBarleyNames(T3ize(rownames(K.gwas2)))

#getting the pheno set up for GWAS function
pheno2 <- pheno23.24.2 %>%
  select(Line, everything())
#factors for BLUEs
pheno2$Location <- factor(pheno2$Location, levels = unique(pheno2$Location))
pheno2$Trial <- factor(pheno2$Trial, levels = unique(pheno2$Trial))
pheno2$Block <- factor(pheno2$Block, levels = unique(pheno2$Block))
pheno2$Line <- factor(pheno2$Line, levels = unique(pheno2$Line))
pheno2$TrialBlock <- factor(paste0(pheno2$Trial, "_", pheno2$Block), levels = unique(paste0(pheno2$Trial, "_", pheno2$Block)))
traits <- c(names(pheno2)[11:18])
#BLUEs
BLUEs2 <- make.BLUEs(pheno2, traits, t = "~TrialBlock + Line")
BLUEs.gwas2 <- BLUEs2
BLUEs.gwas2 <- as.data.frame(BLUEs.gwas2)
BLUEs.gwas2$Line <- rownames(BLUEs.gwas2)
rownames(BLUEs.gwas2) <- NULL
BLUEs.gwas2 <- BLUEs.gwas2[, c("Line", names(BLUEs.gwas2)[1:(ncol(BLUEs.gwas2)-1)])]
#BLUPs
BLUPs2 <- make.BLUPs(pheno2, traits, K.gwas2)
BLUPs.gwas2 <- BLUPs2
BLUPs.gwas2 <- as.data.frame(BLUPs.gwas2)
BLUPs.gwas2$Line <- rownames(BLUPs.gwas2)
rownames(BLUPs.gwas2) <- NULL
BLUPs.gwas2 <- BLUPs.gwas2[, c("Line", names(BLUPs.gwas2)[1:(ncol(BLUPs.gwas2)-1)])]

#getting geno set up for GWAS function
geno.mat2 <- as.matrix(geno)
row.names(geno.mat2) <- cleanBarleyNames(T3ize(row.names(geno.mat2)))
geno.gwas.t2 <- t(geno.mat2 - 1)
lines.keep2 <- intersect(pheno2$Line, colnames(geno.gwas.t2))
geno.gwas.t2 <- geno.gwas.t2[, lines.keep2]
chr2 <- sub("^S(\\d+H)_\\d+$", "\\1", rownames(geno.gwas.t2))
map.pos2 <- as.integer(sub("^S\\d+H_(\\d+)$", "\\1", rownames(geno.gwas.t2)))
geno.gwas2 <- data.frame(Marker = rownames(geno.gwas.t2), Chromosome = chr2, MapPosition = map.pos2, geno.gwas.t2, stringsAsFactors = FALSE, check.names = FALSE) #something weird going on here
colnames(geno.gwas2)[-c(1:3)] <- gsub("\\.", "-", colnames(geno.gwas2)[-c(1:3)])
colnames(geno.gwas2)[-c(1:3)] <- cleanBarleyNames(T3ize(colnames(geno.gwas2)[-c(1:3)]))
all(colnames(K.gwas2) %in% colnames(geno.gwas2))

#GWAS function for GN
BLUPs.gwas.GN2 <- BLUPs.gwas2[, c("Line", "GN")]
reduced <- GWAS(BLUPs.gwas.GN2, geno.gwas2, K = K.gwas2, plot = TRUE)
bonThresh <- -log10(0.05 / length(reduced$Marker)) #5.65
abline(h = bonThresh, col = "firebrick")

sigSNPs <- reduced[reduced$GN >= bonThresh,]
#0 significant SNPs
print(max(reduced$GN))

traits <- BLUPs.gwas2[2:9]
cors <- cor(traits, use = "pairwise.complete.obs")

pdf("../figures/genoCorAll.pdf")
corrplot(cors,
         method = "color",         # color-coded tiles
         type = "lower",           # show only lower triangle
         col = colorRampPalette(c("red", "white", "blue"))(200),  # color scale
         addCoef.col = NA,         # skip showing numbers in boxes
         tl.col = "black",         # text label color
         tl.cex = 0.8,             # label size
         diag = FALSE)  
dev.off()
