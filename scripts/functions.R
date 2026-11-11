#nick's code to clean barley
cleanBarleyNames <- function(lineNames, rmSuffix = FALSE, presoak = FALSE, rmParen = FALSE){
  # lineNames <- c(barley21$Line, barley22$Line)
  suffix <- c("LX", "NA", "LA", "SA", "2R", "6R", "WS", "Dec.", "DEC")
  
  lineNames <- gsub(paste(paste0("\\s*\\((", suffix, ")\\)\\s*$"), collapse = "|"), paste0(" ", paste0("\\", 1:length(suffix), collapse = "")), lineNames)
  if(rmSuffix){
    lineNames <- gsub(paste0("\\s*_*", paste0(suffix, "$"), collapse = "|"), "", lineNames)
  } else {
    lineNames <- gsub(paste0("\\s*", paste0("(", suffix, ")$"), collapse = "|"), paste0("_", paste0("\\", 1:length(suffix), collapse = "")), lineNames)
  }
  if(presoak & rmSuffix) lineNames <- cleanBarleyNames(lineNames, presoak = FALSE)
  lineNames <- gsub("_{2,}", "_", lineNames)
  if(rmParen) lineNames <- gsub("\\(.*\\)$", "", lineNames)
  lineNames
}

#' T3ize function
#'
#' function to reformat line names using T3 rules
#'
#' @param x vector of line names to be reformatted 
#' @return formated line names
#' @details uses T3 rules to reformat line names
#' @examples 
#' T3ize(c("VA96W-403WS*HS", "VA19W-165*", "VA20W-42**", "VA19W-165†", "Jagger", "MCIA Venus", "FFR 555"))
#' @export
T3ize <- function(x){
  rmSpace <- function(xx){
    len <- length(xx)
    if(length(xx) == 1) return(xx)		
    lpos <- nchar(xx[1])
    l <- substr(xx[1], lpos, lpos)
    s <- substr(xx[2], 1, 1)
    bothNum <- grepl("[0-9]", l) & grepl("[0-9]", s)
    bothChar <- grepl("[A-z]", l) & grepl("[A-z]", s)
    if(bothNum | bothChar) adapter <- "_" else adapter <- ""
    xx12 <- paste(xx[1], xx[2], sep = adapter)
    if(len > 2) xx <- c(xx12, xx[3:len]) else xx <- xx12
    rmSpace(xx)
  }
  x <- toupper(x)
  x <- gsub("\"|\'|\`|,", "_", x)
  
  # this might be a bad idea having the dagger match, I just dont know how to recognize that charactre
  endAst <- grep("\\*$|\u2020$", x)
  x[endAst] <- gsub("\\*|\u2020", "", x[endAst])
  midAst <- grep("\\*.+", x)
  x[midAst] <- gsub("\\*", "_", x[midAst])
  x <- gsub("__", "_", x)
  x <- trimws(x)
  xsplit <- strsplit(x, split = " ")
  xsplit <- lapply(xsplit, function(xx) xx[xx != ""])
  sapply(xsplit, rmSpace)
}

#nick's code from CSES 5304

#function to build a QQplot
QQplot <- function(pval, neglog10 = FALSE, ...) {
  E <- rev(-log10(seq(0, 1, length.out = length(pval) + 1)[-1]))
  E <- sort(-log10(runif(length(pval), 0, 1)))
  if(!neglog10) pval <- -log10(pval)
  O <- sort(pval)
  lims <- c(0, max(c(E,O)))
  plot(E, O, xlim = lims, ylim = lims, asp = 1, ...)
  abline(0, 1, col = "dodgerblue", ...)
}

#this function from class will make block diagonal matrices
bdiag <- function(...){
  l <- list(...)
  X <- l[[1]]
  if(length(l) > 1){
    for(i in 2:length(l)){
      X2 <- l[[i]]
      U <- cbind(X, matrix(0, nrow(X), ncol(X2)))
      L <- cbind(matrix(0, nrow(X2), ncol(X)), X2)
      cn <- c(colnames(X), colnames(X2))
      rn <- c(rownames(X), rownames(X2))
      X <- rbind(U, L)
      rownames(X) <- rn
      colnames(X) <- cn
    }
  }
  return(X)
}

#making folds for validation
makeFolds <- function(linenames, k){
  nind <- length(linenames)
  ind <- split(sample(linenames, nind, replace = FALSE), rep(1:k, c(rep(round(nind / k), k - 1), nind - ((k - 1) * (round(nind / k)))))) 
}

#simple function I made to clean names that isn't that great but I am keeping it because it is the first one I ever made
clean.names <- function(names) {
  cleaned <- trimws(tolower(names))
  return(cleaned)
}

#function to make a K matrix
make.K <- function(M, scale_center = 2){
  W <- scale(M, scale = FALSE)
  p <- attributes(W)[["scaled:center"]] / scale_center
  WWt <- tcrossprod(W)
  K <- WWt / (2 * crossprod(p, 1-p)[[1]])
  diag(K) <- diag(K) + 0.01
  return(K)
}

#making a BLUEs function
#example would provide make.BLUEs(barley.pheno, traits, f = "~TrialBlock + Line")
make.BLUEs <- function(pheno.df, traits, t){
  #setting up empty BLUEs
  BLUEs <- matrix(NA, nrow = nlevels(pheno.df$Line), ncol = length(traits), dimnames = list(levels(pheno.df$Line), traits))
  #for loop for calculating
  for(trait in traits){
    # print(trait)
    f <- paste0(trait, t)
    fiti <- lm(as.formula(f), data = pheno.df)
    coefi <- dummy.coef(fiti)
    BLUEs[,trait] <- coefi[["Line"]] + coefi[["(Intercept)"]]
  }
  return(BLUEs)
}

#creating a function for BLUPs
#example setup make.BLUPs(barley.pheno, traits, K)
#uses TrialBlock as the argument, may need to go in and manually edit as necessary
#needs rrBLUP
make.BLUPs <- function(pheno.df, traits, K){
  #empty BLUPs matrix
  BLUPs <- matrix(NA, nrow = nlevels(pheno.df$Line), ncol = length(traits), dimnames = list(levels(pheno.df$Line), traits))
  
  #for loop to run through everything
  for(trait in traits){
    #trait = "PlotWeight.g"
    # print(trait)
    df.trait <- droplevels(pheno.df[!is.na(pheno.df[[trait]]), ])
    #design matrix of fixed effects
    if(length(unique(df.trait$TrialBlock)) == 1){
      X <- model.matrix(~1, data = df.trait)
    } else {
      X <- model.matrix(~TrialBlock, data = df.trait)
    }
    #incidence matrix for random line effects
    Z <- model.matrix(~Line-1, data = df.trait)
    colnames(Z) <- gsub("^Line", "", colnames(Z))# this just cleans up the names of Z
    if(!all(colnames(K) == colnames(Z))) {
      stop("I cant fit the model, the column name of Z dont match the column names of K!")
    }
    gblup.fit <- mixed.solve(y = df.trait[[trait]], X = X, Z = Z, K = K)
    BLUPs[, trait] <- gblup.fit$u
  }
  return(BLUPs)
}

#making function to calculate heritability
#example setup would be make.h2(barley.pheno, traits, K)
#the argument to make X uses the TrialBlock argument
#needs rrBLUP
make.h2 <- function(pheno.df, traits, K) {
  h2 <- NULL
  for (trait in traits) {
    #design matrix of fixed effects
    X <- model.matrix(~ TrialBlock, data = pheno.df)
    colnames(X) <- gsub("^TrialBlock", "", colnames(X))
    
    #incidence matrix for random line effects
    Z <- model.matrix(~Line-1, data = pheno.df)
    colnames(Z) <- gsub("^Line", "", colnames(Z))
    
    #check that Z and K match
    all(colnames(K) == colnames(Z))
    
    gblup.fit <- mixed.solve(y = pheno.df[[trait]], X = X, Z = Z, K = K)
    h2[trait] <- gblup.fit$Vu / (gblup.fit$Vu + gblup.fit$Ve)
  }
  return(h2)
}

#broad-sense heritability
#example setup would be make.H2(barley.pheno, traits)
#the argument to make X uses the TrialBlock argument
#needs rrBLUP
make.H2 <- function(pheno.df, traits) {
  H2 <- NULL
  
  for (trait in traits) {
    
    # design matrix of fixed effects
    X <- model.matrix(~ TrialBlock, data = pheno.df)
    colnames(X) <- gsub("^TrialBlock", "", colnames(X))
    
    # incidence matrix for random line effects
    Z <- model.matrix(~ Line - 1, data = pheno.df)
    colnames(Z) <- gsub("^Line", "", colnames(Z))
    
    # identity matrix for broad-sense heritability
    I <- diag(ncol(Z))
    rownames(I) <- colnames(Z)
    colnames(I) <- colnames(Z)
    
    # fit mixed model
    fit <- mixed.solve(
      y = pheno.df[[trait]],
      X = X,
      Z = Z,
      K = I
    )
    
    # broad-sense heritability
    H2[trait] <- fit$Vu / (fit$Vu + fit$Ve)
  }
  
  return(H2)
}

#function for cross validation and prediction accuracy via folds
#example setup would be fold.cross_pred(barley.pheno, BLUPs, traits, K, folds)
#the argument to make X uses the TrialBlock argument
#should already have BLUEs or BLUPs and K made and folds and n.folds but change n.folds as necessary
#needs rrBLUP
fold.val_pred <- function(pheno.df, traits, K, folds, n.folds = 5) {
  # pheno.df = pheno2;  
  #design matrix of fixed effects
  X <- model.matrix(~ TrialBlock, data = pheno.df)
  colnames(X) <- gsub("^TrialBlock", "", colnames(X))
  
  #incidence matrix for random line effects
  pheno.df$Line <- factor(pheno.df$Line, levels = rownames(K))  # ensure matching order
  
  Z <- model.matrix(~Line-1, data = pheno.df)
  colnames(Z) <- gsub("^Line", "", colnames(Z))
  cor.mat <- matrix(NA, nrow = n.folds, ncol = length(traits), dimnames = list(paste0("k", 1:n.folds), traits))
  mean.cor <- list()
  pred.list.full <- list()
  
  for(trait in traits){
    # trait <- "GN"
    fold.preds <- list()
    fold.obs <- list()
    h2 <- numeric(n.folds)
    
    for(k in 1:n.folds){
      # k = 1
      #define train and test indices
      is.test <- pheno.df$Line %in% folds[[k]] 
      test <- which(is.test)
      train <- which(!is.test)
      #test.lines <- levels(pheno.df$Line)[levels(pheno.df$Line) %in% folds[[k]]]
      
      #subset design matrices
      X.train <- as.matrix(X[train, , drop = FALSE])
      Z.train <- as.matrix(Z[train, , drop = FALSE])
      y.train <- pheno.df[[trait]][train]
      
      #fit model
      ncol(Z.train) == nrow(K)
      all(colnames(Z.train) == rownames(K))
      fit <- mixed.solve(y = y.train, X = X.train, Z = Z.train, K = K)
      h2[k] <- fit$Vu / (fit$Vu + fit$Ve)
      # print(paste0("Trait: ", trait, " ; fold: ", k, " ; h2:", h2))
      #pred.list[[k]] <- train.fit$u[names(train.fit$u) %in% test.lines]
      
      #predict test set
      Z.test <- as.matrix(Z[test, , drop = FALSE])
      y.pred <- Z.test %*% fit$u
      y.obs <- pheno.df[[trait]][test]
      
      #store predictions and observations
      fold.preds[[k]] <- as.numeric(y.pred)
      fold.obs[[k]] <- y.obs
      #per-fold correlation
      cor.mat[k, trait] <- cor(y.obs, y.pred, use = "complete.obs")
    }
    
    #aggregate all folds
    all.preds <- unlist(fold.preds)
    all.obs <- unlist(fold.obs)
    mean.cor[[trait]] <- cor(all.obs, all.preds, use = "complete.obs")
    
    #store full prediction list if desired
    pred.list.full[[trait]] <- list(pred = all.preds, obs = all.obs)
  }  
  return(list(cor.mat = cor.mat,
              mean.cor = mean.cor,
              pred.list = pred.list.full))
}

#function for cross validation and prediction accuracy via train and test input
#example setup would be input.cross_pred(barley.pheno, train.set, test.set, BLUEs, traits, K)
#the argument to make X uses the TrialBlock argument and I have to manually figure out the correct columns
#should already have BLUEs or BLUPs and K made and folds and test and train sets
#needs rrBLUP
input.val_pred <- function(pheno.df, train.df, test.df, traits, K) {
  # Ensure consistent levels and TrialBlock construction
  pheno.df$Line <- factor(pheno.df$Line, levels = rownames(K))
  pheno.df$Trial <- factor(pheno.df$Trial)
  pheno.df$Block <- factor(pheno.df$Block)
  pheno.df$TrialBlock <- factor(paste0(pheno.df$Trial, "_", pheno.df$Block))
  
  train.df$TrialBlock <- factor(paste0(train.df$Trial, "_", train.df$Block),
                                levels = levels(pheno.df$TrialBlock))
  test.df$TrialBlock <- factor(paste0(test.df$Trial, "_", test.df$Block),
                               levels = levels(pheno.df$TrialBlock))
  
  # Build fixed effect matrices
  X.train <- as.matrix(model.matrix(~ TrialBlock, data = train.df))
  colnames(X.train) <- gsub("^TrialBlock", "", colnames(X.train))
  X.test <- as.matrix(model.matrix(~ TrialBlock, data = test.df))
  colnames(X.test) <- gsub("^TrialBlock", "", colnames(X.test))
  
  # Build random effect incidence matrices
  Z.train <- model.matrix(~ Line - 1, data = train.df)
  colnames(Z.train) <- gsub("^Line", "", colnames(Z.train))
  Z.test  <- model.matrix(~ Line - 1, data = test.df)
  colnames(Z.test) <- gsub("^Line", "", colnames(Z.test))
  
  # Align Z matrices to K
  align_Z_to_K <- function(Z, K) {
    all_lines <- rownames(K)
    Z.full <- matrix(0, nrow = nrow(Z), ncol = length(all_lines))
    colnames(Z.full) <- all_lines
    rownames(Z.full) <- rownames(Z)
    match_cols <- intersect(colnames(Z), all_lines)
    Z.full[, match_cols] <- Z[, match_cols, drop = FALSE]
    return(Z.full)
  }
  
  Z.train <- align_Z_to_K(Z.train, K)
  Z.test  <- align_Z_to_K(Z.test, K)
  
  # Final sanity check
  stopifnot(ncol(Z.train) == nrow(K), all(colnames(Z.train) == rownames(K)))
  
  # Output containers
  mean.cor <- numeric(length(traits))
  names(mean.cor) <- traits
  h2 <- numeric(length(traits))
  names(h2) <- traits
  test.pred_list <- list()
  
  # Loop over traits
  for (trait in traits) {
    y.train <- train.df[[trait]]
    y.test  <- test.df[[trait]]
    
    # Fit model
    fit <- rrBLUP::mixed.solve(y = y.train, X = X.train, Z = Z.train, K = K)
    h2[trait] <- fit$Vu / (fit$Vu + fit$Ve)
    
    # Predict
    pred.test <- Z.test %*% fit$u
    pred.test <- as.numeric(pred.test)
    test.pred_list[[trait]] <- pred.test
    
    # Correlation (if test phenotypes exist)
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
  
  return(list(
    test.pred = test.pred_list,
    mean.cor = mean.cor,
    h2 = h2
  ))
}


#libraries needed: dplyr, rrblup, lmer
#need TrialBlock: pheno$TrialBlock <- factor(paste(pheno$Trial, pheno$Block, sep = ":"), levels = unique(paste(pheno$Trial, pheno$Block, sep = ":")))
#traits <- c("GrainYield.bu.ac", "TestWeight.lb.bu", "HeadingDate.Julian", "Height.in", "Lodging.0.9", "PowderyMildew.0.9", "Extract", "BetaGlucan", "AlphaAmylase", "DP", "Protein", "SoverT", "FAN")
#traitDIR <- c("high", "high", "low", "low", "low", "low", "high", "low", "high", "high", "target", "target", "target")
#traitWT <- c(0.20, 0.10, 0.07, 0.07, 0.07, 0.07, 0.12, 0.08, 0.07, 0.05, 0.03, 0.03, 0.03)
#traitTAR optional named vector: c(Protein = 11.5, SoverT = 0.45, FAN = 200)
selectionsWithIndex <- function(pheno, geno, traits, traitDIR, traitWT, traitTAR = NULL, top_n = 20, noGeno = TRUE) {
  
  #filter out traits in pheno not specified
  availTraits <- traits[traits %in% colnames(pheno)]
  traitDIR <- traitDIR[traits %in% availTraits]
  traitWT <- traitWT[traits %in% availTraits]
  
  if (!is.null(traitTAR)) {
    traitTAR <- traitTAR[names(traitTAR) %in% availTraits]
  }
  
  if (length(availTraits) == 0) {
    stop("None of the specified traits are present in phenotype data")
  }
  
  #kinship matrix
  geno <- as.matrix(geno)
  K <- make.K(geno)
  
  #drop lines with no genotypes
  if (noGeno) {
    pheno <- pheno %>% filter(Line %in% rownames(geno))
  }
  
  #get number of observations per line
  lineObs <- pheno %>%
    filter(Line != "" & !is.na(Line)) %>%
    group_by(Line) %>%
    summarise(n_obs = n(), .groups = "drop")
  
  #result list for all trait BLUPs/BLUEs
  blup_list <- list()
  
  for (i in seq_along(availTraits)) {
    trait <- availTraits[i]
    direction <- traitDIR[i]
    weight <- traitWT[i]
    
    #drop rows with NA trait
    ptrait <- pheno %>%
      select(Line, TrialBlock, all_of(trait)) %>%
      filter(!is.na(.data[[trait]]))
    
    # Lines with genotype data
    genoLines <- intersect(unique(ptrait$Line), rownames(K))
    noGenoLines <- setdiff(unique(ptrait$Line), rownames(K))
    
    ### genotyped lines (rrBLUP)
    if (length(genoLines) > 1) {
      # Calculate TrialBlock-adjusted trait means
      pgeno <- ptrait %>% filter(Line %in% genoLines)
      pgeno[[trait]] <- as.numeric(as.character(pgeno[[trait]]))
      model <- lmer(formula(paste(trait, "~ TrialBlock + (1|Line)")),
                    data = pgeno, REML = TRUE)
      ranefs <- ranef(model)$Line
      genoBLUPs <- data.frame(Line = rownames(ranefs), value = ranefs[,1])
    } else {
      genoBLUPs <- data.frame(Line = character(0), value = numeric(0))
    }
    
    ### non genotyped lines (lmer BLUEs)
    if (length(noGenoLines) > 0) {
      pnogeno <- ptrait %>% filter(Line %in% noGenoLines)
      pnogeno[[trait]] <- as.numeric(as.character(pnogeno[[trait]]))
      model <- lm(formula(paste(trait, "~ TrialBlock + Line")),
                  data = pnogeno)
      fixed_ef <- coef(model)
      blueVals <- fixed_ef[grepl("^Line", names(fixed_ef))]
      blueLines <- gsub("Line", "", names(blueVals))
      noGenoBLUEs <- data.frame(Line = blueLines, value = as.numeric(blueVals))
    } else {
      noGenoBLUEs <- data.frame(Line = character(0), value = numeric(0))
    }
    
    ### combine all lines for this trait
    traitBLUPs <- bind_rows(genoBLUPs, noGenoBLUEs)
    colnames(traitBLUPs)[2] <- trait
    blup_list[[trait]] <- traitBLUPs
  }
  
  #merge all trait values into one data frame
  BLUPs <- Reduce(function(x, y) full_join(x, y, by = "Line"), blup_list)
  
  ### standardize and compute index
  for (i in seq_along(traits)) {
    trait <- availTraits[i]
    direction <- traitDIR[i]
    
    # Z-score
    if (!is.null(traitTAR) && trait %in% names(traitTAR)) {
      target <- traitTAR[trait]
      BLUPs[[paste0(trait, "_std")]] <- -abs(BLUPs[[trait]] - target)
    } else {
      z <- scale(BLUPs[[trait]])
      BLUPs[[paste0(trait, "_std")]] <- if (direction == "low") -z else z
    }
  }
  
  #compute weighted index
  BLUPs$index <- rowSums(sapply(seq_along(availTraits), function(i) {
    BLUPs[[paste0(availTraits[i], "_std")]] * traitWT[i]
  }), na.rm = TRUE)
  
  #add number of observations
  BLUPs <- left_join(BLUPs, lineObs, by = "Line")
  
  ### rank and return top lines
  selected <- BLUPs %>%
    arrange(desc(index)) %>%
    slice(1:top_n)
  
  return(selected)
}