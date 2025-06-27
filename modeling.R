# Env
set.seed(1257)

lapply(list(
parallel_processing = "doParallel", 
           modeling = "caret",
                idk = "reshape2", 
                 en = "glmnet", 
               loop = "foreach", 
                 rf = "randomForest",
                svm = "e1071",
          roc_curve = "pROC",
               plot = "ggplot2", 
      color_palette = "viridis"
  ), require, character.only = TRUE
)

# pareto_by_column <- function(valuesVector){
#   scaledVector <- rep(0, length(valuesVector))
#   for (i in 1:length(valuesVector)) {
#     scaledVector[i] <- valuesVector[i] - mean(valuesVector) / sqrt(sd(valuesVector))
#   }
#   return(scaledVector)
# }

customTheme <- theme_minimal() + theme(
    axis.title = element_text(
      family = "Inter", 
      color = "#232323", 
      size = 12
    ),
    axis.title.x = element_text(
      margin = ggplot2::margin(
        t = 13, 
        r = 0, 
        b = 0, 
        l = 0
      )
    ),
    axis.title.y = element_text(
      margin = ggplot2::margin(
        t = 0, 
        r = 13, 
        b = 0, 
        l = 0
      )
    ),
    axis.text = element_text(
      family = "Inter", 
      color = "#232323", 
      size = 8
    )
)


path_lsvm <- "0601_svm/060101_linearSVM/" # linearSVM.R 
path_rsvm <- "0601_svm/060102_radialSVM/" # radialSVM.R 
path_en <- "0602_elastic_net/" # newEN.R
path_rf <- "0603_random_forest/" # newRF.R

# Load and Tidy Measurements Data  
## Create Sample's Data
classes <- cbind.data.frame(
  "Samples" = c(
    paste0("a", c(1:3, 5:15, 17:22)), 
    paste0("t", c(1:3, 5:15, 17:22))
  ), 
  "Diag" = rep(c("NK", "WT"), each = 20, length.out = 40),
  "DiagBin" = rep(c(0L, 1L), each = 20, length.out = 40),

  "Risk" = rep(c("High", "High",
                 "Intermediate", "Intermediate", "Intermediate", 
                 "Intermediate", "Intermediate", 
                 "High", "High", "High", "High", "High", 
                 "Intermediate", "Intermediate", "Intermediate", 
                 "Intermediate", "Intermediate", "Intermediate", 
                 "Intermediate", "Intermediate"), 
                 times = 2),
  "RiskBin" = rep(c(1L, 1L,
                    0L, 0L, 0L, 0L, 0L, 
                    1L, 1L, 1L, 1L, 1L, 
                    0L, 0L, 0L, 0L, 0L, 
                    0L, 0L, 0L), 
                    times = 2)
  )


#########################################################################
merged = FALSE  # TRUE = POS+NEG; FALSE = POS/NEG
mode = TRUE # TRUE = Positive; FALSE = Negative
filtered = FALSE # TRUE = Statistically chosen; FALSE = All processed peaks
#########################################################################
measurements <- NULL
esi <- NULL
significantColNames <- NULL
#########################################################################

if (merged == TRUE) {
  ## Merge Ionization Modes
  dbNeg <- as.data.frame(read.csv("../data/revisedDatabase_neg.csv"))
  dbn <- cbind.data.frame(Condition = as.factor(rep("neg", nrow(dbNeg))), dbNeg)
  names(dbn) <- gsub("neg.", "", names(dbn))
  dbn$Compound <- paste0("neg_", dbn$Compound)

  dbPos <- as.data.frame(read.csv("../data/revisedDatabase_modeling_pos.csv"))
  dbp <- cbind.data.frame(Condition = as.factor(rep("pos", nrow(dbPos))), dbPos)
  names(dbp) <- gsub("pos.", "", names(dbp))
  dbp$Compound <- paste0("pos_", dbp$Compound)

  measurements <- rbind.data.frame(dbn, dbp)
  esi <- "merged/"
  significantColNames <- c(
    "Signif.diag", 
    "Signif.risk", 
    "Signif.risk.nk", 
    "Signif.risk.wt"
  )
} else if (mode == TRUE) {
  # ESI+
  measurements <- as.data.frame(read.csv("../data/revisedDatabase_modeling_pos.csv"))
  esi <- "pos/"
  significantColNames <- c(
    "Signif.pos.diag", 
    "Signif.pos.risk", 
    "Signif.pos.risk.nk", 
    "Signif.pos.risk.wt"
  )

} else {
  # ESI-
  measurements <- as.data.frame(read.csv("../data/revisedDatabase_neg.csv"))
  esi <- "neg/"
  significantColNames <- c(
    "Signif.neg.diag", 
    "Signif.neg.risk", 
    "Signif.neg.risk.nk", 
    "Signif.neg.risk.wt"
  )
}

# ACTUAL ----------------------------------------------------------------|

for (i in 1:length(significantColNames)) {
  comparison <- significantColNames[i]

# Format the numerical matrix to be the base data to all the modeling
baseData <- NULL
if (filtered) {
  baseData <- t(
    measurements[
      which(measurements[ , as.character(comparison)] == 1), 
      c(which(colnames(measurements) == "Names"), 
      which(colnames(measurements) == "a1"):ncol(measurements))
    ]
  )  
} else {
  baseData <- t(
    measurements[ , 
      c(which(colnames(measurements) == "Names"), 
      which(colnames(measurements) == "a1"):ncol(measurements))
    ]
  )
  if (esi == "neg/") {
    esi <- "processed_neg/"
  } else if (esi == "pos/") {
    esi <- "processed_pos/"
  } else if (esi == "merged/") {
    esi <- "processed/"
  }
}

colnames(baseData) <- baseData[1, ]
  
# Normalize by log(1 + x) and Pareto/zScore scaling if needed
baseData <- as.data.frame(
  apply(
    # apply(
      apply(baseData[-1, ], 2, as.numeric),
      2, log1p
      #), 2, pareto_by_column
    )
  )
row.names(baseData) <- classes$Samples

  
# Sample data partition
trainIndex <- sample(
        x = 1:nrow(baseData), 
      size = nrow(baseData) * 0.7, 
  replace = FALSE)
  
  ## Train
  trainSetX <- as.matrix(baseData[trainIndex, ])
  
  if (grepl(".diag", comparison) == TRUE) {
    trainSetY <- as.factor(classes$Diag[trainIndex])
    #print("y: Diag")
  } else {
    trainSetY <- as.factor(classes$Risk[trainIndex])
    #print("y: Risk")
  }  

  ## Test
  testSetX <- as.matrix(baseData[-trainIndex, ])
  
  if (grepl(".diag", comparison) == TRUE) {
    testSetY <- as.factor(classes$Diag[-trainIndex])
    #print("yHat: Diag")
  } else {
    testSetY <- as.factor(classes$Risk[-trainIndex])
    #print("yHat: Risk")
  } 

# Modeling

  #models <- c("")

  models <- c("lsvm", "rsvm", "en", "rf")

  #print.noquote(models)

  print.noquote(paste0("v-------------------- [ ", comparison, " ]"))
  if (sum(grepl("lsvm", models)) == TRUE){
    source(paste0(path_lsvm, "linearSVM.R"))
    print.noquote(paste0("   Linear SVM    [OK]"))
  }
    
  if (sum(grepl("rsvm", models)) == TRUE){
    source(paste0(path_rsvm, "radialSVM.R"))
    print.noquote(paste0("   Radial SVM    [OK]"))
  }

  if (sum(grepl("en", models)) == TRUE){
    source(paste0(path_en, "elasticNets.R"))
    print.noquote(paste0("   Elastic Net   [OK]"))
  }
  
  if (sum(grepl("rf", models)) == TRUE){
    source(paste0(path_rf, "randomForest.R"))
    print.noquote(paste0("   Random Forest [OK]"))
    print.noquote(paste0(""))
  }

}

exportDataMatrix <- NULL
exportDataMatrix <- cbind.data.frame(classes, baseData)
exportDataMatrix$Partition <- "Training"
exportDataMatrix[-trainIndex, "Partition"] <- "Testing"
exportDataMatrix$Partition <- as.factor(exportDataMatrix$Partition)

#write.csv(exportDataMatrix, "neg_partitionData.csv", row.names = FALSE)
#write.csv(exportDataMatrix, "pos_partitionData.csv", row.names = FALSE)
