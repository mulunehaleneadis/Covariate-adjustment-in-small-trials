#Cross-fitted post-LASSO G-computation estimator
## Data generating mechanism
```r
# Required packages
library(MASS)
library(LaplacesDemon)
library(glmnet)
library(arm)

# Number of covariates
p <- 15
# Sample size
n <- 50

# Covariates
generate_cov_matrix <- function(p) {
  outer(1:p, 1:p, function(i, j) 0.1^abs(i - j))
}

Sigma <- generate_cov_matrix(p)
mu <- rep(0, p)
X <- mvrnorm(n, mu = mu, Sigma = Sigma)

# Coefficients
gamma <- (-1)^(1:p) / sqrt(1:p)

# Linear predictor
f_X <- X %*% gamma

# Treatment assignment
A <- rbinom(n, 1, 0.5)

# Outcome
pi <- invlogit(0.25 + 0.75 * A + f_X)
Y <- rbinom(n, 1, pi)

# Data set
data <- data.frame(Y, A, X)
```

## Create cross-fitting folds

```r
fold_ids <- sample(rep(1:5, length.out = dim(data)[1]))
folds <- split(1:dim(data)[1], fold_ids)
```

## Initialize storage for predictions

```r
Q_0_all <- numeric(0)
Q_1_all <- numeric(0)
Y_all <- numeric(0)
A_all <- numeric(0)
```

## LASSO selection, outcome regression, and targeting step

```r
for (i in seq_along(folds)) {
  test_indices <- folds[[i]]
  train_indices <- setdiff(seq_len(nrow(data)), test_indices)

  # Split data into training and test sets
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  X_train <- X[train_indices, , drop = FALSE]
  Y_train <- Y[train_indices]
  A_train <- A[train_indices]

  # Construct the design matrix
  X_train_df <- data.frame(A_train, X_train)
  X_main <- model.matrix(~ ., data = X_train_df)[, -1]
  colnames(X_main) <- make.names(colnames(X_main))

  # Do not penalize treatment
  penalty_factor <- ifelse(colnames(X_main) == "A_train", 0, 1)

  # Fit logistic LASSO model
  lasso_model <- cv.glmnet(X_main,Y_train,alpha = 1,family = "binomial", penalty.factor = penalty_factor)

  # Identify variables selected by LASSO
  LASSO_coef <- as.matrix(coef(lasso_model, s = "lambda.min"))
  nonzero_vars <- rownames(LASSO_coef)[LASSO_coef != 0]
  selected_vars <- nonzero_vars[nonzero_vars %in% colnames(X_main)]

  # Refit the outcome regression (Bayesian GLM) using selected variables
  X_selected <- X_main[, selected_vars, drop = FALSE]
  dta_selected <- data.frame(Y_train, X_selected)
  if (!"A_train" %in% names(dta_selected)) {dta_selected$A_train <- A_train}
  formula_fit <- as.formula(paste("Y_train ~", paste(names(dta_selected)[-1], collapse = "+")))
  model_fit <- bayesglm(formula_fit, data = dta_selected, family = binomial(link = "logit"), x = TRUE)

  # Test data sets under A = 0 and A = 1
  existing_vars <- selected_vars[selected_vars %in% colnames(test_data)]
  if (length(existing_vars) == 0) {
    test_data_0 <- data.frame(A_train = rep(0, nrow(test_data)))
    test_data_1 <- data.frame(A_train = rep(1, nrow(test_data)))
    new_data <- data.frame(A_train = test_data$A)
  } else {
    test_data_0 <- test_data[, existing_vars, drop = FALSE]
    test_data_1 <- test_data[, existing_vars, drop = FALSE]
    test_data_0$A_train <- 0
    test_data_1$A_train <- 1
    new_data <- test_data[, existing_vars, drop = FALSE]
    new_data$A_train <- test_data$A}
  colnames(test_data_0) <- make.names(colnames(test_data_0))
  colnames(test_data_1) <- make.names(colnames(test_data_1))
  colnames(new_data) <- make.names(colnames(new_data))

  # Predict counterfactual outcome probabilities
  pred_0_i <- predict(model_fit, newdata = test_data_0, type = "response")
  pred_1_i <- predict(model_fit, newdata = test_data_1, type = "response")

  # Obtain predictions under the observed treatment
  pred_i <- as.numeric(predict(model_fit, newdata = new_data, type = "response"))

  # Perform the targeting/update step
  epsilon_i <- coef(glm(test_data$Y ~ test_data$A + offset(qlogis(pred_i)), family = binomial()))
  Q0_i <- plogis(qlogis(pred_0_i) + epsilon_i[[1]])
  Q1_i <- plogis(qlogis(pred_1_i) + epsilon_i[[1]] + epsilon_i[[2]])

  # Store fold-specific predictions and observed data
  Q_0_all <- c(Q_0_all, Q0_i)
  Q_1_all <- c(Q_1_all, Q1_i)
  Y_all <- c(Y_all, Y[test_indices])
  A_all <- c(A_all, A[test_indices])
}
```

## Estimate the average treatment effect
```r
G_comp <- mean(Q_1_all - Q_0_all)
```

## Construct the influence function
```r
prob <- mean(A)
inf <- A_all / prob * (Y_all - Q_1_all) + Q_1_all - ((1 - A_all) / (1 - prob) * (Y_all - Q_0_all) +Q_0_all)
```

## Estimate the standard error
```r
n <- nrow(data)
se_G_comp <- sqrt(var(inf) / n)
```

## Construct the 95% confidence interval
```r
CI <- c(lower = G_comp - 1.96 * se_G_comp, upper = G_comp + 1.96 * se_G_comp)
```
