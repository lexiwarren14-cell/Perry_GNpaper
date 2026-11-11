# setwd("")
library(readxl)
library(fieldMapR)
source("functions.R")
gn <- as.data.frame(read_excel("../HartwickGNanalysis95_2023_2R.xlsx"))
gn$Line <- T3ize(cleanBarleyNames(gn$Line))

K <- as.matrix(read.table("../BarleyAddCovWithChecks2023.txt", check.names = FALSE))
dim(K)

gn$Line[!gn$Line %in% colnames(K)]

gn$Line <- factor(gn$Line, levels = unique(gn$Line))

gnData <- gn[!is.na(gn$GN), ]
gnNoData <- gn[is.na(gn$GN), ]

gnData$Line[!gnData$Line %in% colnames(K)]

gnGeno <- droplevels(gn[gn$Line %in% colnames(K), ])
gnGenoNoData <- gnGeno[is.na(gnGeno$GN), ]
gnGenoData <- gnGeno[!is.na(gnGeno$GN), ]

dim(gnGenoData)


length(unique(gnGenoData$Location))

library(asreml)

Kgn <- K[rownames(K) %in% gnGeno$Line, colnames(K) %in% gnGeno$Line]

# fit <- asreml(GN~1, random = ~Line, data = gnGeno)
fit <- asreml(GN~1, random = ~vm(Line, Kgn), data = gnGeno)
summary(fit)$varcomp

pred <- predict(fit, classify = "vm(Line, Kgn)", aliased = TRUE)
head(pred$pvals, 50)

blup <- pred$pvals
names(blup)[1] <- "Line"
blup$predicted <- !blup$Line %in% gnGenoData$Line
blup <- blup[order(blup$predicted.value), ]
blup[c("predicted.value", "std.error")] <- round(blup[c("predicted.value", "std.error")], 3)

write.csv(blup, file = "GNpredictionsSorted.csv", row.names = FALSE)

vc <- summary(fit)$varcomp
sigmaG <- vc[grep("Line", rownames(vc)), "component"]
sigmaG / sum(vc$component)
pevD <- pred$pvals$std.error^2
mean(1 - pevD / sigmaG)

mean(1-blup$std.error[blup$predicted]^2/sigmaG)
mean(1-blup$std.error[!blup$predicted]^2/sigmaG)

blup[blup$predicted, ]

predGN <- blup[blup$predicted, ]

batchSize <- 16
nPerBatch <- batchSize / 2
predGNhi <- predGN[1:batchSize, ]
predGNlo <- predGN[nrow(predGN):{nrow(predGN) - batchSize + 1}, ]

predGNhi$batch <- rep(c(1, 2), each = nPerBatch)
predGNlo$batch <- rep(c(1, 2), each = nPerBatch)

predGNhilo <- rbind(predGNhi, predGNlo)
names(predGNhilo)[1] <- "Line"
newSamplesToRun <- merge(gn, predGNhilo, by = "Line")
newSamplesToRun <- newSamplesToRun[!duplicated(newSamplesToRun$Line), ]

dim(newSamplesToRun)
write.csv(newSamplesToRun, file = "predictedGNsamplesToRun.csv", row.names = FALSE)


