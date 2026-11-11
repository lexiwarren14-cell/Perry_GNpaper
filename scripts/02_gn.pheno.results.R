#libraries and functions
# setwd("")

library(dplyr)
library(rrBLUP)
library(ggplot2)
library(corrplot)

source("functions.R")

#read in files
phenoDataAll <- read.csv("../phenoDataAll.csv", check.names = FALSE)
predicted <- read.csv("../predictedGNsamplesToRun.csv", check.names = FALSE)
predicted$Line <- cleanBarleyNames(T3ize(predicted$Line))
allPred <- read.csv("../GNpredictionsSorted_2R.csv", check.names = FALSE)

#make predicted histogram
pdf("../figures/predHist.pdf")
allPred$Group <- "Middle"
allPred$Group[1:16] <- "Low"
allPred$Group[(nrow(allPred) - 15):nrow(allPred)] <- "High"
allPred$Group <- factor(allPred$Group, levels = c("Low", "Middle", "High"))
dist <- ggplot(allPred, aes(x = reorder(Line, predicted.value), y = predicted.value, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Low" = "blue", "Middle" = "grey", "High" = "red")) +
  labs(x = "Line",
       y = "Predicted GN Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(labels = allPred$LineLabel)
print(dist)
dev.off()

#subset by year
results <- phenoDataAll[phenoDataAll$Year == 2023, ]
gn24.all <- phenoDataAll[phenoDataAll$Year == 2024, ]

predicted$headType <- "2"

#combining predicted and results datasets
common.cols <- intersect(names(predicted), names(results))
combine <- merge(results, predicted, by = common.cols, all.x = TRUE)

combine.all <- combine[!is.na(combine$GN), ]
combine.all <- combine.all[order(combine.all$GN), ]

lines <- c("DH162310", "DH170472", "VA20MFHB-18DH349", "VA21M-24", "CALYPSO")
hist.data <- combine.all[combine.all$Line %in% lines, ]
hist.data$GN <- as.numeric(hist.data$GN)
hist.data$Line <- c("GN-0 Line 1", "GN-0 Line 2", "Low-GN Line 1", "Low-GN Line 2", "Calypso")

pdf("../figures/lineHist.pdf")
ggplot(hist.data, aes(x = Line, y = GN, fill = Line)) +
  geom_bar(stat = "identity") +
  labs(x = "Line", y = "GN Value (ppm)", pch = 16) +
  theme_minimal() +
  theme(legend.position = "none", axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))
dev.off()

pdf("../figures/dist23.pdf")
hist(combine.all$GN, xlab = "GN Value (ppm)", main = "Distribution of GN in 2023")
dev.off()
print(max(combine.all$GN)) #5.3
print(min(combine.all$GN)) #0.3
print(mean(combine.all$GN)) #2.5

#filtering out the batches and getting in order
combine.filter <- combine[!is.na(combine$batch), ]
combine.filter$GN <- as.numeric(combine.filter$GN)
combine.filter <- combine.filter %>% 
  arrange(batch, predicted.value)

#assigning low to first 8 of each batch and high to last 8
rows.per.batch <- nrow(combine.filter) / 4
combine.filter$high.low <- c(rep("low", rows.per.batch), rep("high", rows.per.batch))

#simple average of the GN value of each batch
averages <- aggregate(GN ~ high.low, data = combine.filter, FUN = mean)
print(averages)
#high = 2.5875
#low = 1.7375

#min a max for each batch
low.min <- min(combine.filter[combine.filter$high.low == "low", ]$GN)
low.max <- max(combine.filter[combine.filter$high.low == "low", ]$GN)
high.min <- min(combine.filter[combine.filter$high.low == "high", ]$GN)
high.max <- max(combine.filter[combine.filter$high.low == "high", ]$GN)
print(low.min) #1.0
print(low.max) #3.2
print(high.min) #0.8
print(high.max) #4.1

#simple t-test
t_test <- t.test(GN ~ high.low, data = combine.filter)
print(t_test)
#t-value = 3.2843
#df = 25.172
#p-value = 0.003
#95% confidence interval: 0.317 to 1.383
combine.filter$high.low <- factor(combine.filter$high.low, levels = c("low", "high"))

#boxplot
pdf("../figures/boxplotGN.jpg")
boxplot(GN ~ high.low, data = combine.filter,
        ylab = "GN Value",
        xlab = "GN Bi-Directional Genomic Selection",
        ylim = c(0, 5),
        main = "Boxplot of GN Values",
        col = c("blue", "red"),
        names = c("Low", "High"))
segments(1, 4.5, 1, 4.8)
segments(1, 4.8, 1.45, 4.8)
segments(1.55, 4.8, 2, 4.8)
segments(2, 4.5, 2, 4.8)
text(1.5, 4.8, labels = "**")
text(1.5, 4.5, labels = "p = 0.003")
dev.off()

#want to know gn-zero, low-producers and high producers
zeros <- combine.all[combine.all$GN < 0.5, ]
low <- combine.all[combine.all$GN >= 0.5 & combine.all$GN <= 1.5, ]
high <- combine.all[combine.all$GN > 1.5, ]

ultra <- combine.all[combine.all$GN <= 0.5, ]

#plot predictions against actual values
pdf("../figures/predVSmeas.pdf")
plot(combine.filter$predicted.value, combine.filter$GN,
     xlab = "Predicted", ylab = "Measured",
     pch = 16,
     col = ifelse(combine.filter$high.low == "low", "blue", "red"))
dev.off()
cor(combine.filter$predicted.value, combine.filter$GN)

combineLow <- combine.filter[combine.filter$high.low == "low", ]
combineHigh <- combine.filter[combine.filter$high.low == "high", ]
plot(combineLow$predicted.value, combineLow$GN, xlab = "Predicted", ylab = "Actual", main = "GN Levels (ppm) of the Low Set", pch = 16)
plot(combineHigh$predicted.value, combineHigh$GN, xlab = "Predicted", ylab = "Actual", main = "GN Levels (ppm) of the High Set", pch = 16)

#2024 phenotype analyses
gn24 <- gn24.all

pheno24 <- gn24[!is.na(gn24$GN), ]
pheno24 <- gn24[!is.na(gn24$moisture), ]
pheno24 <- gn24[!is.na(gn24$fineExtract), ]
pheno24 <- gn24[!is.na(gn24$SP), ]
pheno24 <- gn24[!is.na(gn24$FAN), ]
pheno24 <- gn24[!is.na(gn24$BG), ]
pheno24 <- gn24[!is.na(gn24$AA), ]
pheno24 <- gn24[!is.na(gn24$DP), ]
pheno24$Line <- cleanBarleyNames(T3ize(pheno24$Line))

pheno24$GN <- as.numeric(pheno24$GN)
monResults <- pheno24[!is.na(pheno24$GN), ]
#combining hartwick 24 and montana 25 datasets
combine25 <- merge(results, monResults, by = "Line", all.x = TRUE)
combine25$GN.x <- as.numeric(combine25$GN.x)

combine.all25 <- combine25[!is.na(combine25$GN.x) | !is.na(combine25$GN.y), ]
combine.all25.GN <- combine.all25[, c("Line", "GN.x", "GN.y", "plot_name.x", "plot_name.y")]
combine.all25.GN <- combine.all25.GN[order(combine.all25.GN$GN.x), ]

twoYrGN <- combine.all25.GN[!is.na(combine.all25.GN$GN.x) & !is.na(combine.all25.GN$GN.y), ]
unique(twoYrGN$Line)
cor(twoYrGN$GN.x, twoYrGN$GN.y) #0.73
pdf("../figures/twoYrCor.pdf")
plot(twoYrGN$GN.x, twoYrGN$GN.y, xlab = "2023 GN", ylab = "2024 GN", main = "Phenotypic Correlation", pch = 16)
text(0.8, 3.9, labels = "r = 0.73")
dev.off()

cor(twoYrGN$GN.x, twoYrGN$GN.y)

pdf("../figures/dist24.pdf")
hist(monResults$GN, xlab = "GN Value (ppm)", main = "Distribution of GN in 2024")
dev.off()

#trait correlations across 2023 and 2024
gn24 <- gn24.all
gn24$GN <- as.numeric(gn24$GN)
gn24 <- gn24[!is.na(gn24$GN), ]
gn24 <- gn24[!is.na(gn24$moisture), ]
gn24 <- gn24[!is.na(gn24$fineExtract), ]
gn24 <- gn24[!is.na(gn24$SP), ]
gn24 <- gn24[!is.na(gn24$FAN), ]
gn24 <- gn24[!is.na(gn24$BG), ]
gn24 <- gn24[!is.na(gn24$AA), ]
gn24 <- gn24[!is.na(gn24$DP), ]
gn24$Line <- cleanBarleyNames(T3ize(gn24$Line))

gn23 <- results
gn23 <- gn23[!is.na(gn23$GN), ]
gn23$Line <- cleanBarleyNames(T3ize(gn23$Line))

#combining 23 and 24 gn data
common.cols <- intersect(names(gn23), names(gn24))
pheno23.24 <- rbind(gn23[, common.cols], gn24[, common.cols])

traits <- pheno23.24[11:18]
cors <- cor(traits, use = "pairwise.complete.obs")
pdf("../figures/phenoCorAll.pdf")
corrplot(cors,
         method = "color",         # color-coded tiles
         type = "lower",           # show only lower triangle
         col = colorRampPalette(c("red", "white", "blue"))(200),  # color scale
         addCoef.col = NA,         # skip showing numbers in boxes
         tl.col = "black",         # text label color
         tl.cex = 0.8,             # label size
         diag = FALSE)  
dev.off()

#correlation matrix
cors <- cor(traits, use = "pairwise.complete.obs")

#create empty matrix for p-values
p.mat <- matrix(NA, ncol(traits), ncol(traits))
colnames(p.mat) <- colnames(traits)
rownames(p.mat) <- colnames(traits)

#calculate pairwise correlation p-values
for(i in 1:ncol(traits)){
  for(j in 1:ncol(traits)){
    
    test <- cor.test(traits[[i]],
                     traits[[j]],
                     method = "pearson")
    
    p.mat[i,j] <- test$p.value
  }
}

#view p-values
round(p.mat, 4)
p.mat["GN", "SP"]
p.mat["GN", "FAN"]
p.mat["GN", "DP"]
p.mat["GN", "BG"]

pdf("../figures/bgVSgn2324.pdf")
plot(pheno23.24$BG, pheno23.24$GN, xlab = "Beta-Glucan", ylab = "Glycosidic Nitrile", pch = 16)
dev.off()
cor(pheno23.24$BG, pheno23.24$GN)
#-0.21

feed23 <- gn23[grepl("^VA\\d+B-", gn23$Line) | gn23$Line == "SECRETARIAT", ] #should be empty just validating for errors
malt23 <- gn23[!(gn23$Line %in% feed23$Line), ]

pdf("../figures/bgVSgn23.pdf")
plot(gn23$BG, gn23$GN, xlab = "Beta-Glucan", ylab = "Glycosidic Nitrile", pch = 16, main = "2023")
dev.off()

pdf("../figures/bgVSgn24.pdf")
plot(gn24$BG, gn24$GN, xlab = "Beta-Glucan", ylab = "Glycosidic Nitrile", pch = 16, main = "2024")
dev.off()