
set.seed(7121)
# lapply(list(
#              modeling = "glmnet", 
#                   idk = "reshape2", 
#                  loop = "foreach", 
#   parallel_processing = "doParallel", 
#                  plot = "ggplot2", 
#         color_palette = "viridis"), 
#   require, character.only = TRUE)

# Reset Variables
alphaValues <- NULL
coreCluster <- NULL
search <- NULL
cv <- NULL
cvProperties <- NULL
bestAlpha <- NULL
bestLambda <- NULL
enModel <- NULL
pred <- NULL
pred_prob <- NULL
coefs <- NULL
coefDF <- NULL
coefPlot <- NULL
infLimit <- NULL
supLimit <- NULL

# Elastic Net Model (Ridge: alpha = 0; Elastic: 0 < alpha < 1; LASSO: alpha = 1)
# Manual-ish Cross-Validation to search for the optimal Lambda value
alphaValues <- seq(0.1, 0.9, 0.05)
# It loops (foreach; for loop) in parallel processing (registerDoParallel) a 5-fold CV, 
#   with each alpha value of seq(...) to find the best lambda
coreCluster <- makeCluster(detectCores() - 1)
registerDoParallel(coreCluster)

search <- foreach(i = alphaValues, .combine = rbind ) %dopar% { 
  cv <- glmnet::cv.glmnet(trainSetX, trainSetY,
    family = "binomial", alpha = i, nfold = 5, 
    type.measure = "deviance", paralle = TRUE
  )
  data.frame(
    cvError = cv$cvm[cv$lambda == cv$lambda.1se], 
    lambda.1se = cv$lambda.1se, alpha = i
  )
} 

stopCluster(coreCluster)

cvProperties <- search[search$cvError == min(search$cvError), ] #// find the minimum cross-validated mean error (best model)
bestAlpha <- cvProperties$alpha
bestLambda <- cvProperties$lambda.1se

write.csv(search, paste0(path_en, esi, comparison,"_enCV.csv"), row.names = FALSE)

# Fit the Final Model With the Best Hyperparameters
enModel <- glmnet(
  trainSetX, trainSetY,
  alpha = bestAlpha, 
  lambda = bestLambda,
  family = "binomial", 
  probability = TRUE
)

## Predict testSetX based on the trained model enModel
pred <- predict(enModel, testSetX, s = bestLambda, type = "class", probability = TRUE)
pred_prob <- predict(enModel, testSetX, s = bestLambda, type = "response", probability = TRUE)
write.csv(cbind.data.frame(
  Samples = rownames(pred),
  PredictedClass = as.factor(pred), 
  Probability = as.numeric(pred_prob)
), 
  paste0(path_en, esi, comparison,  "_enPredictions.csv"), row.names = FALSE)

# Plot for paper
## Variable Importance Coefficient (continuous data)
coefs <- NULL
coefs <- unlist(as.matrix(coef(enModel, s = bestLambda)))
#coefs <- cbind(Features = rownames(coefs), coefs)
#coefs$features[grep("^X", coefs$features)] <- gsub("X", "", coefs$features[grep("^X", coefs$features)])
#coefs$features <- gsub("m.z", "m/z", coefs$features)

coefDF <- cbind.data.frame(Features = as.character(rownames(coefs)),
                        Coefficient = as.numeric(unlist(coefs)))

coefDF <- coefDF[-1,]

pattern <- "\\d+\\.\\d+_\\d+\\.\\d+(\\w+/\\w+|n)"
coefDF$Identification <- 1
coefDF$Identification[grep(pattern, coefDF$Features)] <- 0

write.csv(coefDF, paste0(path_en, esi, comparison, "_enCoefficients.csv"), row.names = FALSE)                        

coefPlot <- coefDF[which(coefDF$Coefficient < -0.001 | coefDF$Coefficient > 0), ] # remove Coef == 0
limits <- range(coefPlot$Coefficient)

## Identified significative compounds
png(paste0(path_en, esi, comparison, "_enCoefficientsIdentified.png"), 2300, 3200, res = 300) # Signif.neg.diag
print.noquote(ggplot(coefPlot[coefPlot$Identification == 1, ], 
  aes(x = reorder(Features, Coefficient), 
      y = Coefficient, fill = Coefficient)) +
  geom_bar(stat = "identity", position = "dodge", na.rm = TRUE) +      
  geom_hline(yintercept = 0, lwd = 0.2, color = "#232323") +
  scale_y_continuous(limits = c(limits[1] - (0.2 * limits[1]), limits[2] + (0.2 * limits[2]))) + 
  coord_flip() + # Editable: limits = c(-0.2, 0.5)
  scale_fill_viridis() +
  geom_text(aes(label = format(round(Coefficient, 4), scientific = FALSE), 
                hjust = ifelse(Coefficient < 0, 1.2, -0.2)), 
                size.unit = "pt", size = 8, color = "#77767b") +
  labs(
    # title = "Elastic Net Classification", 
    # subtitle = "Variable Importance (Coefficients)", 
          x = "Predictors", 
          y = "Coefficient Value") +
  theme_minimal() + theme(element_text(family = "Inter"),
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.1),
    title = element_text(color = "#232323", size = 12, hjust = 0.5),
    axis.title = element_text(color = "#232323", size = 10),
    axis.text = element_text(color = "#77767b"),
    legend.key.size = unit(0.6, units = "cm"),
    legend.title = element_text(color = "#232323"),
    legend.text = element_text(color = "#77767b")))
dev.off()

## All significative compounds
png(paste0(path_en, esi, comparison, "_enCoefficients.png"), 2300, 3200, res = 300) # Signif.neg.diag
print.noquote(ggplot(coefPlot, 
  aes(x = reorder(Features, Coefficient), 
      y = Coefficient, fill = Coefficient)) +
  geom_bar(stat = "identity", position = "dodge", na.rm = TRUE, width = 0.5) +      
  geom_hline(yintercept = 0, lwd = 0.2, color = "#232323") +
  scale_y_continuous(limits = c(limits[1] - 0.3, limits[2] + 0.5)) + coord_flip() + # Editable: limits = c(-0.2, 0.5)
  scale_fill_viridis() +
  geom_text(aes(label = format(round(Coefficient, 4), scientific = FALSE), 
                hjust = ifelse(Coefficient < 0, 1.2, -0.2)), 
                size.unit = "pt", size = 8, color = "#77767b") +
  labs(
    # title = "Elastic Net Classification", 
    # subtitle = "Variable Importance (Coefficients)", 
          x = "Predictors", 
          y = "Coefficient Value") +
  theme_minimal() + theme(element_text(family = "Inter"),
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.1),
    title = element_text(color = "#232323", size = 12, hjust = 0.5),
    axis.title = element_text(color = "#232323", size = 10),
    axis.text = element_text(color = "#77767b"),
    legend.key.size = unit(0.6, units = "cm"),
    legend.title = element_text(color = "#232323"),
    legend.text = element_text(color = "#77767b")))
dev.off()

## Density plot significative compounds
png(paste0(path_en, esi, comparison, "_enCoefficientsDensity.png"), 1800, 2300, res = 300) # Signif.neg.diag
ggplot(coefPlot, aes(x = Coefficient)) + 
geom_histogram(fill = "#E3E6FFFF", binwidth = 0.00025) + 
geom_density(kernel = "gaussian", color = "#404788", fill = "#404788", alpha = 0.25) +
labs(y = "Density") + 
#coord_cartesian(xlim = c(0, max(coefPlot$Coefficient) + 0.0025)) + 
customTheme
dev.off()

## CV Plot for Ridge/LASSO
#plot(search$cvError, search$lambda.1se)

# ROC Curve Evaluation
source("roc.R")
png(paste0(path_en, esi, comparison, "_enROC.png"), width = 2000, height = 2000, res = 300)
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

write.csv(rocData, paste0(path_en, esi, comparison, "_enROC.csv"))
write.csv(confMatrix_table, paste0(path_en, esi, comparison, "_enConfusionMatrix.csv"))
write.csv(confMatrix_metrics, paste0(path_en, esi, comparison, "_enMetrics.csv"))
