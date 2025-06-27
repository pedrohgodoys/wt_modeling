
set.seed(1257)
# lapply(list(
#       roc = "pROC", 
#      plot = "ggplot2", 
#   palette = "viridis"), 
#   require, character.only = TRUE
# )

rocCurve <- NULL
rocData <- NULL
confMatrix <- NULL
TP <- NULL
TN <- NULL
FP <- NULL
FN <- NULL
modelMetrics <- NULL
catThreshold <- 0.7 # tirei do cu

# ROC Curve 
if (grepl(".diag", comparison) == TRUE) {
    rocCurve <- roc(
       response = as.factor(testSetY), 
      predictor = as.numeric(pred_prob[, 1]), 
         levels = c("NK", "WT"), 
      direction = ">"
    )
  } else {
    rocCurve <- roc(
       response = as.factor(testSetY), 
      predictor = as.numeric(pred_prob[, 1]), 
         levels = c("Intermediate", "High"), 
      direction = ">"
    )
  } 

rocData <- NULL
  rocData <- data.frame(
    Specificity = rev(rocCurve$specificities), 
    Sensitivity = rev(rocCurve$sensitivities)
  )

# Confusion matrix and model parameters calculation
confMatrix <- confusionMatrix(as.factor(pred), as.factor(testSetY))
confMatrix_table <- confMatrix$table
confMatrix_metrics <- as.data.frame(t(confMatrix$byClass))

# confMatrix <- table(Predicted = pred, Actual = testSetY)

# if (dim(confMatrix)[1] == 1){
#   confMatrix <- rbind(confMatrix, Missed = 0) 
  
# }

#   ## Values
#   TP <- as.numeric(confMatrix[2,2])
#   TN <- as.numeric(confMatrix[1,1])
#   FP <- as.numeric(confMatrix[2,1])
#   FN <- as.numeric(confMatrix[1,2])
  
#   ## Model Metrics
#   modelMetrics <- cbind.data.frame(
#     Accuracy = (TP + TN) / sum(confMatrix),
#     Sensitivity = TP / (TP + FN),
#     Specificity = TN / (TN + FP),
#     Precision = TP / (TP + FP),
#     F1score = 2 * ((TP / (TP + FP)) * (TP / (TP + FN))) / 
#       ((TP / (TP + FP)) + (TP / (TP + FN)))
#   )
#   modelMetrics[is.nan(unlist(modelMetrics))] <- 0
