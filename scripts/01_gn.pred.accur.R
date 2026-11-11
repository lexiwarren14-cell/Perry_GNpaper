#libraries and functions
# setwd("")

library(rrBLUP)
library(dplyr)

source("functions.R")


#read in files
geno <- read.table("../BarleyMarkersConsensusWithChecks2024.txt")
phenoDataAll <- read.csv("../phenoDataAll.csv", check.names = FALSE)

#kinship matrix
K.all <- make.K(geno)
colnames(K.all) <- cleanBarleyNames(T3ize(colnames(K.all)))
rownames(K.all) <- cleanBarleyNames(T3ize(rownames(K.all)))

#subset by year
gn23.all <- phenoDataAll[phenoDataAll$Year == 2023, ]
gn24.all <- phenoDataAll[phenoDataAll$Year == 2024, ]


#2023 analysis
pheno <- gn23.all[!is.na(gn23.all$GN), ]

#match phenotype lines to marker lines
K <- K.all
all(colnames(K) %in% pheno$Line)
all(pheno$Line %in% colnames(K))

no.geno <- unique(pheno$Line[!pheno$Line %in% colnames(K)])
pheno2 <- pheno[!pheno$Line %in% no.geno, ]

pheno2$Line <- factor(pheno2$Line, levels = unique(pheno2$Line))
lines <- levels(pheno2$Line)
K <- K[lines, lines]
all(lines == colnames(K))

#factors
pheno2$Location <- factor(pheno2$Location, levels = unique(pheno2$Location))
pheno2$Trial <- factor(pheno2$Trial, levels = unique(pheno2$Trial))
pheno2$Block <- factor(pheno2$Block, levels = unique(pheno2$Block))
pheno2$Line <- factor(pheno2$Line, levels = unique(pheno2$Line))
pheno2$Year <- factor(pheno2$Year, levels = unique(pheno2$Year))
pheno2$TrialBlock <- factor(
  paste0(pheno2$Trial, "_", pheno2$Block),
  levels = unique(paste0(pheno2$Trial, "_", pheno2$Block))
)

traits <- c(names(pheno2)[11:18])

#BLUEs
BLUEs23 <- make.BLUEs(pheno2, traits, t = "~TrialBlock + Line")

#BLUPs
BLUPs23 <- matrix(
  NA,
  nrow = nlevels(pheno2$Line),
  ncol = length(traits),
  dimnames = list(levels(pheno2$Line), traits)
)

for (trait in traits) {
  df.trait <- droplevels(pheno2[!is.na(pheno2[[trait]]), ])
  
  if (nrow(df.trait) > 0) {
    if (length(unique(df.trait$TrialBlock)) == 1) {
      X <- model.matrix(~1, data = pheno2)
    } else {
      X <- model.matrix(~TrialBlock, data = pheno2)
    }
    
    Z <- model.matrix(~Line - 1, data = pheno2)
    colnames(Z) <- gsub("^Line", "", colnames(Z))
    
    if (!all(colnames(K) == colnames(Z))) {
      stop("I cant fit the model, the column name of Z dont match the column names of K!")
    }
    
    gblup.fit <- mixed.solve(
      y = pheno2[[trait]],
      X = X,
      Z = Z,
      K = K
    )
    
    BLUPs23[, trait] <- gblup.fit$u
  } else {
    message(paste("Trait", trait, "has insufficient data for BLUP calculation"))
  }
}

# Heritability.
h2.2023 <- make.h2(pheno2, traits, K)
print(h2.2023)
#GN = 0.25

#10-fold cross-validation
n.folds <- 10
n.reps <- 10
set.seed(0)
mean.cor.list <- vector("list", n.reps)

for (i in seq_len(n.reps)) {
  folds <- makeFolds(levels(pheno2$Line), n.folds)
  gn.pred <- fold.val_pred(pheno2, traits, K, folds)
  mean.cor.list[[i]] <- as.numeric(gn.pred$mean.cor)
  
  if (!is.null(names(gn.pred$mean.cor))) {
    names(mean.cor.list[[i]]) <- names(gn.pred$mean.cor)
  }
  
  cat("Completed repetition", i, "of", n.reps, "\n")
}

mean.cor.mat.2023 <- do.call(rbind, mean.cor.list)
avg.mean.cor.2023 <- colMeans(mean.cor.mat.2023, na.rm = TRUE)
print(mean.cor.mat.2023)
print(avg.mean.cor.2023)
#GN = 0.36


#2024 analysis
pheno24 <- gn24.all[!is.na(gn24.all$GN), ]
pheno24$GN <- as.numeric(pheno24$GN)
pheno24 <- pheno24[!is.na(pheno24$GN), ]

K <- K.all
all(colnames(K) %in% pheno24$Line)
all(pheno24$Line %in% colnames(K))

no.geno <- unique(pheno24$Line[!pheno24$Line %in% colnames(K)])
pheno24.2 <- pheno24[!pheno24$Line %in% no.geno, ]

pheno24.2$Line <- factor(pheno24.2$Line, levels = unique(pheno24.2$Line))
lines24 <- levels(pheno24.2$Line)
K <- K[lines24, lines24]
all(lines24 == colnames(K))

pheno24.2$Location <- factor(pheno24.2$Location, levels = unique(pheno24.2$Location))
pheno24.2$Trial <- factor(pheno24.2$Trial, levels = unique(pheno24.2$Trial))
pheno24.2$Block <- factor(pheno24.2$Block, levels = unique(pheno24.2$Block))
pheno24.2$Line <- factor(pheno24.2$Line, levels = unique(pheno24.2$Line))
pheno24.2$TrialBlock <- factor(
  paste0(pheno24.2$Trial, "_", pheno24.2$Block),
  levels = unique(paste0(pheno24.2$Trial, "_", pheno24.2$Block))
)

traits <- "GN"

BLUEs24 <- make.BLUEs(pheno24.2, traits, t = "~TrialBlock + Line")

h2.2024 <- make.h2(pheno24.2, traits, K)
print(h2.2024)
#GN = 0.37

n.folds <- 10
n.reps <- 10
set.seed(0)
mean.cor.list <- vector("list", n.reps)

for (i in seq_len(n.reps)) {
  folds <- makeFolds(levels(pheno24.2$Line), n.folds)
  gn.pred <- fold.val_pred(pheno24.2, traits, K, folds)
  mean.cor.list[[i]] <- as.numeric(gn.pred$mean.cor)
  
  if (!is.null(names(gn.pred$mean.cor))) {
    names(mean.cor.list[[i]]) <- names(gn.pred$mean.cor)
  }
  
  cat("Completed repetition", i, "of", n.reps, "\n")
}

mean.cor.mat.2024 <- do.call(rbind, mean.cor.list)
avg.mean.cor.2024 <- colMeans(mean.cor.mat.2024, na.rm = TRUE)
print(mean.cor.mat.2024)
print(avg.mean.cor.2024)
#GN = 0.46

BLUPs24 <- matrix(
  NA,
  nrow = nlevels(pheno24.2$Line),
  ncol = length(traits),
  dimnames = list(levels(pheno24.2$Line), traits)
)

for (trait in traits) {
  df.trait <- droplevels(pheno24.2[!is.na(pheno24.2[[trait]]), ])
  
  if (nrow(df.trait) > 0) {
    if (length(unique(df.trait$TrialBlock)) == 1) {
      X <- model.matrix(~1, data = pheno24.2)
    } else {
      X <- model.matrix(~TrialBlock, data = pheno24.2)
    }
    
    Z <- model.matrix(~Line - 1, data = pheno24.2)
    colnames(Z) <- gsub("^Line", "", colnames(Z))
    
    if (!all(colnames(K) == colnames(Z))) {
      stop("I cant fit the model, the column name of Z dont match the column names of K!")
    }
    
    gblup.fit <- mixed.solve(
      y = pheno24.2[[trait]],
      X = X,
      Z = Z,
      K = K
    )
    
    BLUPs24[, trait] <- gblup.fit$u
  } else {
    message(paste("Trait", trait, "has insufficient data for BLUP calculation"))
  }
}


#genotypic correlation between 2023 and 2024 BLUPs
BLUPs23.df <- as.data.frame(BLUPs23)
BLUPs23.df$Line <- row.names(BLUPs23.df)

BLUPs24.df <- as.data.frame(BLUPs24)
BLUPs24.df$Line <- row.names(BLUPs24.df)

BLUPsTwoYr <- merge(BLUPs23.df, BLUPs24.df, by = "Line")
BLUPsTwoYr <- BLUPsTwoYr[BLUPsTwoYr$Line %in% pheno2$Line, ]

cor(BLUPsTwoYr$GN.x, BLUPsTwoYr$GN.y) #0.79

pdf("../figures/genoCorr.pdf")
plot(
  BLUPsTwoYr$GN.x,
  BLUPsTwoYr$GN.y,
  xlab = "2023 GN",
  ylab = "2024 GN",
  main = "Genotypic Correlation",
  pch = 16
)
text(-0.3, 1.0, labels = "r = 0.79")
dev.off()

#combined 2023 and 2024 GN analysis

pheno23.24 <- phenoDataAll[
  phenoDataAll$Year %in% c(2023, 2024) & !is.na(phenoDataAll$GN),
]

K <- K.all
all(colnames(K) %in% pheno23.24$Line)
all(pheno23.24$Line %in% colnames(K))

no.geno <- unique(pheno23.24$Line[!pheno23.24$Line %in% colnames(K)])
pheno23.24.2 <- pheno23.24[!pheno23.24$Line %in% no.geno, ]

pheno23.24.2$Line <- factor(
  pheno23.24.2$Line,
  levels = unique(pheno23.24.2$Line)
)
lines23.24 <- levels(pheno23.24.2$Line)
K <- K[lines23.24, lines23.24]
all(lines23.24 == colnames(K))

pheno23.24.2$Location <- factor(
  pheno23.24.2$Location,
  levels = unique(pheno23.24.2$Location)
)
pheno23.24.2$Trial <- factor(
  pheno23.24.2$Trial,
  levels = unique(pheno23.24.2$Trial)
)
pheno23.24.2$Block <- factor(
  pheno23.24.2$Block,
  levels = unique(pheno23.24.2$Block)
)
pheno23.24.2$Line <- factor(
  pheno23.24.2$Line,
  levels = unique(pheno23.24.2$Line)
)
pheno23.24.2$Year <- factor(
  pheno23.24.2$Year,
  levels = unique(pheno23.24.2$Year)
)
pheno23.24.2$TrialBlock <- factor(
  paste0(pheno23.24.2$Trial, "_", pheno23.24.2$Block),
  levels = unique(paste0(pheno23.24.2$Trial, "_", pheno23.24.2$Block))
)

traits <- "GN"

BLUEs23.24 <- make.BLUEs(
  pheno23.24.2,
  traits,
  t = "~TrialBlock + Line + Year"
)

h2.23.24 <- make.h2(pheno23.24.2, traits, K)
print(h2.23.24)
#GN = 0.43

BLUPs23.24 <- matrix(
  NA,
  nrow = nlevels(pheno23.24.2$Line),
  ncol = length(traits),
  dimnames = list(levels(pheno23.24.2$Line), traits)
)

for (trait in traits) {
  df.trait <- droplevels(pheno23.24.2[!is.na(pheno23.24.2[[trait]]), ])
  
  if (nrow(df.trait) > 0) {
    if (length(unique(df.trait$TrialBlock)) == 1) {
      X <- model.matrix(~1, data = pheno23.24.2)
    } else {
      X <- model.matrix(~TrialBlock, data = pheno23.24.2)
    }
    
    Z <- model.matrix(~Line - 1, data = pheno23.24.2)
    colnames(Z) <- gsub("^Line", "", colnames(Z))
    
    if (!all(colnames(K) == colnames(Z))) {
      stop("I cant fit the model, the column name of Z dont match the column names of K!")
    }
    
    gblup.fit <- mixed.solve(
      y = pheno23.24.2[[trait]],
      X = X,
      Z = Z,
      K = K
    )
    
    BLUPs23.24[, trait] <- gblup.fit$u
  } else {
    message(paste("Trait", trait, "has insufficient data for BLUP calculation"))
  }
}

n.folds <- 10
n.reps <- 10
set.seed(0)
mean.cor.list <- vector("list", n.reps)

for (i in seq_len(n.reps)) {
  folds <- makeFolds(levels(pheno23.24.2$Line), n.folds)
  gn.pred <- fold.val_pred(pheno23.24.2, traits, K, folds)
  mean.cor.list[[i]] <- as.numeric(gn.pred$mean.cor)
  
  if (!is.null(names(gn.pred$mean.cor))) {
    names(mean.cor.list[[i]]) <- names(gn.pred$mean.cor)
  }
  
  cat("Completed repetition", i, "of", n.reps, "\n")
}

mean.cor.mat.23.24 <- do.call(rbind, mean.cor.list)
avg.mean.cor.23.24 <- colMeans(mean.cor.mat.23.24, na.rm = TRUE)
print(mean.cor.mat.23.24)
print(avg.mean.cor.23.24)
#GN = 0.36

#combined-year malting-quality analysis
#retain only records with complete values for the original eight traits.
traits <- c("GN", "moisture", "fineExtract", "SP", "FAN", "BG", "AA", "DP")
traits <- traits[traits %in% names(phenoDataAll)]

pheno2 <- phenoDataAll
for (trait in traits) {
  pheno2 <- pheno2[!is.na(pheno2[[trait]]), ]
}

K <- K.all
no.geno <- unique(pheno2$Line[!pheno2$Line %in% colnames(K)])
pheno2 <- pheno2[!pheno2$Line %in% no.geno, ]

pheno2$Line <- factor(pheno2$Line, levels = unique(pheno2$Line))
lines <- levels(pheno2$Line)
K <- K[lines, lines]
all(lines == colnames(K))

pheno2$Location <- factor(pheno2$Location, levels = unique(pheno2$Location))
pheno2$Trial <- factor(pheno2$Trial, levels = unique(pheno2$Trial))
pheno2$Block <- factor(pheno2$Block, levels = unique(pheno2$Block))
pheno2$Line <- factor(pheno2$Line, levels = unique(pheno2$Line))
pheno2$Year <- factor(pheno2$Year, levels = unique(pheno2$Year))
pheno2$TrialBlock <- factor(
  paste0(pheno2$Trial, "_", pheno2$Block),
  levels = unique(paste0(pheno2$Trial, "_", pheno2$Block))
)

BLUEsCombined <- make.BLUEs(pheno2, traits, t = "~TrialBlock:Year + Line")

h2.combined <- make.h2(pheno2, traits, K)
print(h2.combined)
#GN = 0.42

n.folds <- 10
n.reps <- 10
set.seed(0)
mean.cor.list <- vector("list", n.reps)

for (i in seq_len(n.reps)) {
  folds <- makeFolds(levels(pheno2$Line), n.folds)
  gn.pred <- fold.val_pred(pheno2, traits, K, folds)
  mean.cor.list[[i]] <- as.numeric(gn.pred$mean.cor)
  
  if (!is.null(names(gn.pred$mean.cor))) {
    names(mean.cor.list[[i]]) <- names(gn.pred$mean.cor)
  }
  
  cat("Completed repetition", i, "of", n.reps, "\n")
}

mean.cor.mat.combined <- do.call(rbind, mean.cor.list)
avg.mean.cor.combined <- colMeans(mean.cor.mat.combined, na.rm = TRUE)
print(mean.cor.mat.combined)
print(avg.mean.cor.combined)
#GN = 0.36

BLUPsCombined <- matrix(
  NA,
  nrow = nlevels(pheno2$Line),
  ncol = length(traits),
  dimnames = list(levels(pheno2$Line), traits)
)

for (trait in traits) {
  df.trait <- droplevels(pheno2[!is.na(pheno2[[trait]]), ])
  
  if (nrow(df.trait) > 0) {
    if (length(unique(df.trait$TrialBlock)) == 1) {
      X <- model.matrix(~1, data = pheno2)
    } else {
      X <- model.matrix(~TrialBlock, data = pheno2)
    }
    
    Z <- model.matrix(~Line - 1, data = pheno2)
    colnames(Z) <- gsub("^Line", "", colnames(Z))
    
    if (!all(colnames(K) == colnames(Z))) {
      stop("I cant fit the model, the column name of Z dont match the column names of K!")
    }
    
    gblup.fit <- mixed.solve(
      y = pheno2[[trait]],
      X = X,
      Z = Z,
      K = K
    )
    
    BLUPsCombined[, trait] <- gblup.fit$u
  } else {
    message(paste("Trait", trait, "has insufficient data for BLUP calculation"))
  }
}


#train on 2023 and predict 2024
train.df <- droplevels(pheno2[pheno2$Year == 2023, ])
test.df <- droplevels(pheno2[pheno2$Year == 2024, ])

train.df$Line <- factor(train.df$Line, levels = rownames(K))
test.df$Line <- factor(test.df$Line, levels = rownames(K))

train.df$TrialBlock <- factor(paste0(train.df$Trial, "_", train.df$Block))
test.df$TrialBlock <- factor(paste0(test.df$Trial, "_", test.df$Block))

#fixed-effect matrix used for model fitting.
X.train <- model.matrix(~TrialBlock, data = train.df)

#random-effect incidence matrices.
Z.train <- model.matrix(~Line - 1, data = train.df)
colnames(Z.train) <- gsub("^Line", "", colnames(Z.train))

Z.test <- model.matrix(~Line - 1, data = test.df)
colnames(Z.test) <- gsub("^Line", "", colnames(Z.test))

#align Z.train to K
Z.train.full <- matrix(
  0,
  nrow = nrow(Z.train),
  ncol = length(rownames(K)),
  dimnames = list(rownames(Z.train), rownames(K))
)
match.cols <- intersect(colnames(Z.train), rownames(K))
Z.train.full[, match.cols] <- Z.train[, match.cols, drop = FALSE]
Z.train <- Z.train.full

#align Z.test to K
Z.test.full <- matrix(
  0,
  nrow = nrow(Z.test),
  ncol = length(rownames(K)),
  dimnames = list(rownames(Z.test), rownames(K))
)
match.cols <- intersect(colnames(Z.test), rownames(K))
Z.test.full[, match.cols] <- Z.test[, match.cols, drop = FALSE]
Z.test <- Z.test.full

stopifnot(all(colnames(Z.train) == rownames(K)))
stopifnot(all(colnames(Z.test) == rownames(K)))

mean.cor <- numeric(length(traits))
names(mean.cor) <- traits
h2.2023.predict.2024 <- numeric(length(traits))
names(h2.2023.predict.2024) <- traits
test.pred.list <- list()

for (trait in traits) {
  y.train <- train.df[[trait]]
  y.test <- test.df[[trait]]
  
  fit <- rrBLUP::mixed.solve(
    y = y.train,
    X = X.train,
    Z = Z.train,
    K = K
  )
  
  h2.2023.predict.2024[trait] <- fit$Vu / (fit$Vu + fit$Ve)
  
  pred.test <- as.numeric(Z.test %*% fit$u)
  test.pred.list[[trait]] <- pred.test
  
  if (all(is.na(y.test))) {
    mean.cor[trait] <- NA
    message("No phenotypic data in test set for trait: ", trait)
  } else {
    valid <- complete.cases(y.test, pred.test)
    
    if (sum(valid) > 1) {
      mean.cor[trait] <- cor(y.test[valid], pred.test[valid])
    } else {
      mean.cor[trait] <- NA
      message("Not enough complete cases to calculate correlation for trait: ", trait)
    }
  }
}

print(mean.cor) #GN = 0.46
print(h2.2023.predict.2024) #GN = 0.25


#GN prediction plots
#store GN predictions in the 2024 test data
test.df$pred_GN <- test.pred.list$GN

test.df$pred_GN <- test.df$pred_GN + 0.5533706

#separate lines observed in both years from lines observed only in 2024
unObs24 <- test.df[!test.df$Line %in% train.df$Line, ]
Obs24 <- test.df[test.df$Line %in% train.df$Line, ]

cor(unObs24$GN, unObs24$pred_GN, use = "complete.obs") #GN = 0.36
cor(Obs24$GN, Obs24$pred_GN, use = "complete.obs") #GN = 0.73

pdf('../figures/unObs24.pdf')
plot(
  unObs24$pred_GN,
  unObs24$GN,
  xlab = "Predicted",
  ylab = "Measured",
  main = "Measured vs Predicted GN Levels in 2024",
  pch = 16
)
abline(lm(unObs24$GN ~ unObs24$pred_GN))
text(-0.7, 3.5, labels = "r = 0.36")
dev.off()

pdf("../figures/obs24.pdf")
plot(
  Obs24$pred_GN,
  Obs24$GN,
  xlab = "Predicted",
  ylab = "Measured",
  main = "Measured vs Predicted GN Levels in 2024",
  pch = 16
)
abline(lm(Obs24$GN ~ Obs24$pred_GN))
text(-0.7, 3.5, labels = "r = 0.73")
dev.off()

pdf("../figures/measVSpred24.pdf")
xlim <- range(0, 1.5)
ylim <- range(0, 5.5)

plot(
  test.df$pred_GN,
  test.df$GN,
  type = "n",
  xlab = "Predicted GN",
  ylab = "Measured GN",
  main = "Measured vs Predicted GN Levels in 2024",
  xlim = xlim,
  ylim = ylim
)

points(Obs24$pred_GN, Obs24$GN, pch = 16, col = "blue")
points(unObs24$pred_GN, unObs24$GN, pch = 16, col = "red")

legend(
  "topleft",
  legend = c(
    "Observed in 2023 and 2024 (r = 0.73)",
    "Observed only in 2024 (r = 0.36)"
  ),
  col = c("blue", "red"),
  pch = 16,
  lwd = 2,
  bty = "n"
)
dev.off()
