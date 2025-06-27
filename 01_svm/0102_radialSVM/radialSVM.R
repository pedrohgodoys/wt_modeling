
set.seed(7354)
# lapply(list(
# modeling = "caret",
#      svm = "e1071",
# plotting = "ggplot2", "doParallel"), 
# require, character.only = TRUE
# )

# Reset Variable
coreCluster <- NULL
svmTrained <- NULL
svmModel <- NULL
pred <- NULL
pred_prob <- NULL


# 5-fold Cross-validation and Model
coreCluster <- makeCluster(detectCores() - 1) # Leave one core free
registerDoParallel(coreCluster)
svmTrained <- caret::train(x = trainSetX, y = trainSetY,
    method = "svmRadial", 
  trControl = trainControl(
    method = "cv", 
    number = 5
  )
)
stopCluster(coreCluster)

write.csv(as.data.frame(svmTrained$results), paste0(path_rsvm, esi, comparison, "_rsvmCV.csv"), row.names = FALSE)

# Fit the actual SVM model
svmModel <- svm(
  x = trainSetX, y = trainSetY, 
  kernel = "radial", 
  probability = TRUE
)
#write.csv(as.data.frame(svmModel), paste0(exportAddr, "_rsvmModel.csv"), row.names = FALSE)

# Predict the testSet
## Classifications
pred <- predict(
  svmModel, testSetX, 
  probability = FALSE, type = "class"
)

## Probabilities
pred_prob <- NULL
pred_prob <- as.data.frame(
  attr(
    predict(
      svmModel, 
      testSetX, 
      probability = TRUE
    ), 
    "probabilities"
  )[, 2]
)

write.csv(
  cbind.data.frame(
    Samples = classes$Samples[-trainIndex],
    PredictedProb = pred_prob[, 1], 
    predictedClass = pred
  ), 
  paste0(path_rsvm, esi, comparison, "_rsvmPredictions.csv"), 
  row.names = FALSE
)

# ROC Curve Evaluation
source("roc.R")
png(paste0(path_rsvm, esi, comparison, "_rsvmROC.png"), width = 2000, height = 2000, res = 300)
print(
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

write.csv(rocData, paste0(path_rsvm, esi, comparison, "_rsvmROC.csv"))
write.csv(confMatrix_table, paste0(path_rsvm, esi, comparison, "_rsvmConfusionMatrix.csv"))
write.csv(confMatrix_metrics, paste0(path_rsvm, esi, comparison, "_rsvmMetrics.csv"))
