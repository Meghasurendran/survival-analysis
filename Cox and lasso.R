library(survminer)
library(survival)
library(glmnet)
library(pROC)
library(MASS)
data=read.csv(file.choose())
length(data)
head(data)
colnames(data)
data$Tumour=(toupper(data$Tumour))
km_fit=survfit(Surv(Time_in_month,Status) ~ 1,data=data)
options(max.print = 10000)
summary(km_fit)
plot(km_fit,xlab = "Time in days",ylab = "Survival probability",main="Kaplan Meier Survival Curve")

surv_object = Surv(time = data$Time_in_month, event = data$Status)
#
# Variables to test
variables <- c("Age_Category","Sex","PST","Histology","Grade","Tumour","Nodal_Status","Metastasis_status","Stage","Treatment_Given_Prior_to_Registration_at_RI","Treatment_Given","Tobacco_chewing","smoking_habit","alcohol_consumption") # Replace with your variable names
logrank_results <- data.frame(Variable = character(),
                              p_value = numeric(),
                              stringsAsFactors = FALSE)
for (var in variables) {
  formula <- as.formula(paste("surv_object ~", var))
  tryCatch({
    log_rank_test <- survdiff(formula, data = data)
    p_value <- 1 - pchisq(log_rank_test$chisq, df = length(log_rank_test$n) - 1)
    logrank_results <- rbind(logrank_results, data.frame(Variable = var, p_value = p_value))
  }, error = function(e) {
    logrank_results <- rbind(logrank_results, data.frame(Variable = var, p_value = NA))
  })
}
print(logrank_results)

#cox ph

# List of covariates
covariates <- c("Age_Category","Sex","PST","Histology","Grade","Tumour","Nodal_Status","Metastasis_status","Stage","Treatment_Given_Prior_to_Registration_at_RI","Treatment_Given","Tobacco_chewing","smoking_habit","alcohol_consumption")

# Initialize a list to store results
univariate_results <- list()

# Loop through each covariate and fit a univariate Cox model
for (var in covariates) {
  formula <- as.formula(paste("Surv(Time_in_month, Status) ~", var))
  cox_model <- coxph(formula, data = data)
  univariate_results[[var]] <- summary(cox_model)
}

# Print results for all variables
for (var in names(univariate_results)) {
  cat("\n--- Univariate Cox Model for:", var, "---\n")
  print(univariate_results[[var]])
  # Create the formula for the multivariate model
  formula <- as.formula(paste("Surv(Time_in_month, Status) ~", paste(var, collapse = " + ")))
  # Fit the multivariate Cox model
  cox_model <- coxph(formula, data = data)
  ph_test <- cox.zph(cox_model)
  # Print PH Test Results
  print(ph_test)
  # Extract p-values for the covariate and global test
  individual_p <- round(ph_test$table[1, 3], 4)  # Column 3 is the p-value
  global_p <- round(ph_test$table[nrow(ph_test$table), 3], 4)
  # Plot Schoenfeld Residuals with p-value annotation
  plot(ph_test)
  
  # Annotate plot with p-values
  mtext(side = 3, line = 0.5, adj = 0, 
        text = paste("Individual p-value:", individual_p, 
                     "| Global p-value:", global_p), 
        cex = 0.8)
}

# Define combinations of variables (start with the full set of variables)
variable_combinations <- list(
  "Sex + PST + Histology + Grade + Tumour + Nodal_Status + Metastasis_status + Stage +Treatment_Given_Prior_to_Registration_at_RI + Treatment_Given + Tobacco_chewing + smoking_habit + alcohol_consumption" = 
    "Sex + PST + Histology + Grade + Tumour + Nodal_Status + Metastasis_status + Stage +Treatment_Given_Prior_to_Registration_at_RI+ Treatment_Given + Tobacco_chewing + smoking_habit + alcohol_consumption"
)

#Stepwise backward elimination

# Initialize an empty data frame to store the AIC results at each step
aic_results <- data.frame(
  Step = integer(),
  Variables = character(),
  AIC = numeric()
)

# Loop through the combinations and calculate AIC for each model
for (combo_name in names(variable_combinations)) {
  # Create the formula for the model using as.formula
  formula <- as.formula(paste("Surv(Time_in_month, Status) ~", variable_combinations[[combo_name]]))
  
  # Fit the initial Cox model using the full set of variables
  model <- coxph(formula, data = data)
  
  # Perform backward stepwise selection and store AIC at each step
  step_model <- stepAIC(model, direction = "backward", trace = TRUE)
  
  # Capture and print the AIC at each step
  for (i in 1:length(step_model$history$AIC)) {
    cat("Step", i, "AIC:", step_model$history$AIC[i], "Variables:", step_model$history$terms[i], "\n")
  }
}

# List of covariates
covariates <- c("PST", "Grade", "Tumour","Nodal_Status","Treatment_Given", "smoking_habit")
# Create the formula for the multivariate model
formula <- as.formula(paste("Surv(Time_in_month, Status) ~", paste(covariates, collapse = " + ")))
# Fit the multivariate Cox model
cox_model_multivariate_backward <- coxph(formula, data = data)
# Summary of the multivariate Cox model
summary(cox_model_multivariate_backward)

X <- model.matrix(~ Age_Category+Sex+PST+Histology+Grade+Tumour+ Nodal_Status+Metastasis_status+Stage+Treatment_Given_Prior_to_Registration_at_RI+ Treatment_Given+Tobacco_chewing+ smoking_habit +alcohol_consumption- 1, data = data)
y <- Surv(data$Time_in_month, data$Status)
lasso_cv <- cv.glmnet(X, y, family = "cox", alpha = 1)
cat("Optimal lambda:", lasso_cv$lambda.min, "\n")
# Extract coefficients at optimal lambda
lasso_coefs <- coef(lasso_cv, s = lasso_cv$lambda.min)  
print(lasso_coefs)
plot(lasso_cv, main = "Cross-validation for Lasso Cox Model")

# Select columns that match partial names for variables of interest
selected_vars <- grep("Age_Category|PST|Grade|Tumour|Nodal_Status|Metastasis_status|Stage|Treatment_Given|Tobacco_chewing|smoking_habit|alcohol_consumption", 
                      colnames(X), value = TRUE)
final_X <- X[, selected_vars]

# Fit final model with selected variables
final_lasso_model <- glmnet(final_X, y, family = "cox", alpha = 1, lambda = lasso_cv$lambda.min)
# Compute risk scores for the final model
final_lasso_risk_scores <- predict(final_lasso_model, newx = final_X, s = lasso_cv$lambda.min, type = "link")

# Fit a Cox model using risk scores
final_lasso_cox_model <- coxph(Surv(data$Time_in_month, data$Status) ~ final_lasso_risk_scores)
final_lasso_cindex <- summary(final_lasso_cox_model)$concordance[1]
cat("C-index for Final Lasso Cox Model:", final_lasso_cindex, "\n")
# Generate ROC curve and compute AUC
roc_final_lasso <- roc(data$Status, final_lasso_risk_scores)
final_lasso_auc <- auc(roc_final_lasso)
cat("AUC for Final Lasso Cox Model:", final_lasso_auc, "\n")

# C-index for the Cox model
cox_cindex <- summary(cox_model_multivariate_backward)$concordance[1]
cat("C-index for Cox model:", cox_cindex, "\n")
# Extracting the linear predictors from the Cox model
cox_pred <- cox_model_multivariate_backward$linear.predictors
# Generate the ROC curve for the Cox model
roc_cox <- roc(data$Status, cox_pred)
# Calculate AUC for the Cox model
auc_cox <- auc(roc_cox)
cat("AUC for Cox model:", auc_cox, "\n")


