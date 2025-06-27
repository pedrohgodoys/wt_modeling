
set.seed(1346)
# lapply(list(
#   "randomForest", "caret", 
#   "ggplot2", "viridis"), 
#   require, character.only = TRUE
# )

# Reset Variables
coreCluster <- NULL
rfModel <- NULL
pred <- NULL
pred_prob <- NULL
rfModelImp <- NULL

n_accur <- function(vector) {
  naccur <- as.data.frame(table(vector))
  
  return(naccur[naccur$Freq > 1, ])
}

# Training Model
coreCluster <- makeCluster(detectCores() - 1)
registerDoParallel(coreCluster)
rfModel <- caret::train(x = trainSetX, y = trainSetY, method = "rf", trControl = trainControl("cv", 5))
stopCluster(coreCluster)

# Predictions
pred <- predict(rfModel, testSetX, type = "raw")
pred_prob <- predict(rfModel, testSetX, type = "prob")

write.csv(cbind.data.frame(
  Samples = classes$Samples[-trainIndex],
  PredictedClass = as.factor(pred), 
  Probability = as.numeric(pred_prob[, 1])), 
  paste0(path_rf, esi, comparison,"_rfPredictions.csv"), row.names = FALSE)

# Importance Plot
rfModelImp <- randomForest(x = trainSetX, y = trainSetY, importance = TRUE)
varImpPlot(rfModelImp)

n_accur(rownames(rfModelImp$importance))

## Organize data - Duplicates values, CAUTION!
rfImportance <- NULL
rfImportance <- cbind.data.frame(
  Compound = rownames(rfModelImp$importance), 
  rfModelImp$importance)
rfImportance$Identification <- ifelse(grepl(pattern, rfImportance[,1]), 0, 1)

n_accur(rownames(rfImportance))

write.csv(
  rfImportance, 
  paste0(path_rf, esi, comparison, "_rfImportance.csv"), 
  row.names = FALSE
)

pattern <- "^X\\d+\\.\\d+_\\d+\\.\\d+(\\w+/\\w+|n)"

## Density distribution of all compounds
png(paste0(path_rf, esi, comparison, "_rfCoefficientDensity.png"), width = 2000, height = 1200, res = 300)
print.noquote(
  ggplot(rfImportance[rfImportance$MeanDecreaseAccuracy > 0,], aes(x = MeanDecreaseAccuracy)) + 
    geom_histogram(fill = "#E3E6FFFF", binwidth = 0.00025) + 
    geom_density(kernel = "gaussian", color = "#404788", fill = "#404788", alpha = 0.25) +
    labs(y = "Density") + 
      coord_cartesian(xlim = c(0, max(rfImportance$MeanDecreaseAccuracy) + 0.0025)) + 
        customTheme)
      dev.off()

## Identified compounds
rfImpKnown <- rfImportance[which(rfImportance$MeanDecreaseAccuracy > 0 & rfImportance$Identification == 1), ]


n_accur(rfImpKnown$Compound)


png(paste0(path_rf, esi, comparison, "_rfCoefficientIdentified.png"), width = 2300, height = 4000, res = 300)
print.noquote(
  ggplot(
    rfImpKnown, 
    aes(x = MeanDecreaseAccuracy,
        y = reorder(Compound, 
                    MeanDecreaseAccuracy))
  ) +
    geom_bar(
      stat = "identity", 
      position = "dodge", 
      fill = "#404788", 
      width = 0.4) +
    xlim(
      0, 
      max(rfImpKnown$MeanDecreaseAccuracy) + (0.2 * max(rfImpKnown$MeanDecreaseAccuracy))
        ) +
    labs(y = "Compounds") +
    geom_text(
      aes(
        label = format(round(MeanDecreaseAccuracy, 4), scientific = FALSE), 
        hjust = ifelse(MeanDecreaseAccuracy < 0.025, -0.1, 1.2)
  ), 
    size.unit = "pt", size = 10, color = "black") +
   customTheme
  
)
dev.off()

# ROC Curve Evaluation
source("roc.R")
png(paste0(path_rf, esi, comparison, "_rfROC.png"), width = 2000, height = 2000, res = 300)
print.noquote(
  ggplot(rocData, 
    aes(x = 1 - Specificity, 
        y = Sensitivity)) +
  geom_area(alpha = 0.1, fill = "#404788") +
  geom_abline(slope = 1, linetype = "dashed", color = "#AAAAAAFF") +
  geom_line(linewidth = 1.5, color = "#404788") + 
  geom_point(size = 3, color = "#404788") +
  annotate("text", x = 0.9, y = 0.15, 
  fontface = "bold", label = paste("AUC =", round(auc(rocCurve), 4)), size = 5, color = "#404788") +
  customTheme
)
dev.off()

write.csv(rocData, paste0(path_rf, esi, comparison, "_rfROC.csv"))
write.csv(confMatrix_table, paste0(path_rf, esi, comparison, "_rfConfusionMatrix.csv"))
write.csv(confMatrix_metrics, paste0(path_rf, esi, comparison, "_rfMetrics.csv"))
