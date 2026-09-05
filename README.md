# Cross-Fitted LASSO Outcome Regression with Targeting for ATE Estimation

This code implements an estimation procedure for the **Average Treatment Effect (ATE)** using:

* 5-fold cross-fitting
* LASSO variable selection
* Bayesian logistic outcome regression
* A targeting/update step
* Influence-function-based standard errors
* A 95% confidence interval

## Required Packages

```r
library(glmnet)
library(arm)
```

If needed, install the packages first:

```r
install.packages(c("glmnet", "arm"))
```

## Required Data Objects

The code assumes that the following objects are already available in the R environment:

* `data`: data frame containing the observed data
* `X`: matrix or data frame containing baseline covariates
* `Y`: binary outcome variable
* `A`: binary treatment indicator coded as `0` or `1`

The number of rows in `data` and `X`, and the lengths of `Y` and `A`, should be equal.

For reproducibility, it is recommended to set a random seed before creating the folds:

```r
set.seed(123)
```

## Step 1: Create Cross-Fitting Folds

The observations are randomly divided into five folds.

```r
fold_ids <- sample(rep(1:5, length.out = dim(data)[1]))
folds <- split(1:dim(data)[1], fold_ids)
```

## Step 2: Initialize Storage for Cross-Fitted Predictions

```r
Q_0_all <- numeric(0)
Q_1_all <- numeric(0)
Y_all <- numeric(0)
A_all <- numeric(0)
```

## Step 3: LASSO Selection, Outcome Regression, and Targeting

For each fold, the remaining observations are used as the training sample and the held-out fold is used as the validation sample.

LASSO is first used for variable selection. The treatment variable is not penalized.

The selected variables are then used to refit the outcome regression using `bayesglm`.

Predictions are obtained under treatment values `A = 0` and `A = 1`, followed by a targeting/update step.

```r
for (i in seq_along(folds)) {

  test_indices <- folds[[i]]
  train_indices <- setdiff(seq_len(nrow(data)), test_indices)

  # Split data into training and validation sets
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]

  X_train <- X[train_indices, , drop = FALSE]
  Y_train <- Y[train_indices]
  A_train <- A[train_indices]

  # Construct the design matrix
  X_train_df <- data.frame(A_train, X_train)

  X_main <- model.matrix(
    ~ .,
    data = X_train_df
  )[, -1]

  colnames(X_main) <- make.names(colnames(X_main))

  # Do not penalize treatment
  penalty_factor <- ifelse(
    colnames(X_main) == "A_train",
    0,
    1
  )

  # Fit the logistic LASSO model
  lasso_model <- cv.glmnet(
    X_main,
    Y_train,
    alpha = 1,
    family = "binomial",
    penalty.factor = penalty_factor
  )

  # Identify variables selected by LASSO
  LASSO_coef <- as.matrix(
    coef(
      lasso_model,
      s = "lambda.min"
    )
  )

  selected_vars <- rownames(LASSO_coef)[LASSO_coef != 0]

  selected_vars <- selected_vars[
    selected_vars %in% colnames(X_main)
  ]

  # Refit the outcome regression using selected variables
  X_selected <- X_main[
    ,
    selected_vars,
    drop = FALSE
  ]

  dta_selected <- data.frame(
    Y_train,
    X_selected
  )

  # Ensure treatment remains in the model
  if (!"A_train" %in% names(dta_selected)) {
    dta_selected$A_train <- A_train
  }

  formula_fit <- as.formula(
    paste(
      "Y_train ~",
      paste(
        names(dta_selected)[-1],
        collapse = "+"
      )
    )
  )

  model_fit <- bayesglm(
    formula_fit,
    data = dta_selected,
    family = binomial(link = "logit"),
    x = TRUE
  )

  # Create test data under A = 0 and A = 1
  existing_vars <- selected_vars[
    selected_vars %in% colnames(test_data)
  ]

  if (length(existing_vars) == 0) {

    test_data_0 <- data.frame(
      A_train = rep(
        0,
        nrow(test_data)
      )
    )

    test_data_1 <- data.frame(
      A_train = rep(
        1,
        nrow(test_data)
      )
    )

    new_data <- data.frame(
      A_train = test_data$A
    )

  } else {

    test_data_0 <- test_data[
      ,
      existing_vars,
      drop = FALSE
    ]

    test_data_1 <- test_data[
      ,
      existing_vars,
      drop = FALSE
    ]

    test_data_0$A_train <- 0
    test_data_1$A_train <- 1

    new_data <- test_data[
      ,
      existing_vars,
      drop = FALSE
    ]

    new_data$A_train <- test_data$A
  }

  colnames(test_data_0) <- make.names(
    colnames(test_data_0)
  )

  colnames(test_data_1) <- make.names(
    colnames(test_data_1)
  )

  colnames(new_data) <- make.names(
    colnames(new_data)
  )

  # Predict counterfactual outcome probabilities
  pred_0_i <- predict(
    model_fit,
    newdata = test_data_0,
    type = "response"
  )

  pred_1_i <- predict(
    model_fit,
    newdata = test_data_1,
    type = "response"
  )

  # Predictions under the observed treatment
  pred_i <- as.numeric(
    predict(
      model_fit,
      newdata = new_data,
      type = "response"
    )
  )

  # Targeting/update step
  epsilon_i <- coef(
    glm(
      test_data$Y ~
        test_data$A +
        offset(qlogis(pred_i)),
      family = binomial
    )
  )

  Q01_i <- plogis(
    qlogis(pred_0_i) +
      epsilon_i[[1]]
  )

  Q11_i <- plogis(
    qlogis(pred_1_i) +
      epsilon_i[[1]] +
      epsilon_i[[2]]
  )

  # Store fold-specific predictions
  Q_0_all <- c(
    Q_0_all,
    Q01_i
  )

  Q_1_all <- c(
    Q_1_all,
    Q11_i
  )

  Y_all <- c(
    Y_all,
    Y[test_indices]
  )

  A_all <- c(
    A_all,
    A[test_indices]
  )
}
```

## Step 4: Estimate the Average Treatment Effect

The ATE is estimated as the average difference between the updated predicted outcomes under treatment and control:

```r
G_comp <- mean(
  Q_1_all - Q_0_all
)
```

Here:

* `Q_1_all` contains updated predicted outcomes under `A = 1`
* `Q_0_all` contains updated predicted outcomes under `A = 0`

## Step 5: Construct the Influence Function

First estimate the marginal probability of treatment:

```r
prob <- mean(A)
```

Then calculate the influence-function values:

```r
inf <- A_all / prob *
  (Y_all - Q_1_all) +
  Q_1_all -
  (
    (1 - A_all) / (1 - prob) *
      (Y_all - Q_0_all) +
      Q_0_all
  )
```

## Step 6: Estimate the Standard Error

```r
n <- nrow(data)

se_G_comp <- sqrt(
  var(inf) / n
)
```

## Step 7: Construct the 95% Confidence Interval

```r
CI <- c(
  lower = G_comp - 1.96 * se_G_comp,
  upper = G_comp + 1.96 * se_G_comp
)
```

## Results

The main quantities of interest are:

```r
G_comp
se_G_comp
CI
```

where:

* `G_comp` is the estimated Average Treatment Effect
* `se_G_comp` is the estimated standard error
* `CI` contains the lower and upper limits of the 95% confidence interval

The results can also be printed as:

```r
cat(
  "ATE:",
  G_comp,
  "\n"
)

cat(
  "Standard Error:",
  se_G_comp,
  "\n"
)

cat(
  "95% Confidence Interval:",
  CI["lower"],
  "to",
  CI["upper"],
  "\n"
)
```

## Notes

* The outcome `Y` is assumed to be binary.
* The treatment `A` is assumed to be binary and coded as `0` or `1`.
* The treatment variable is assigned a penalty factor of `0` in the LASSO model, meaning that treatment is not penalized.
* Five-fold cross-fitting is used so that predictions for each observation are generated from a model that was not trained on that observation.
* The targeting step updates the initial outcome-regression predictions before calculation of the treatment effect.
* The influence function is used to estimate the standard error of the ATE estimator.
* The 95% confidence interval uses the normal approximation:

```text
ATE ± 1.96 × Standard Error
```

* Because `qlogis()` is applied to predicted probabilities, probabilities equal to exactly `0` or `1` may produce infinite values. In practical applications, it may be useful to bound predicted probabilities slightly away from `0` and `1`.
