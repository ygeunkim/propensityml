test_that(
  "Fitting GLM for propensity",
  {
    fit <-
      chemical %>%
      ps_glm(poisox ~ age + sex, data = .)
    expect_s3_class(fit, "propmod")
    expect_s3_class(fit$model, "glm")
    expect_length(fit, 3)
  }
)

test_that(
  "Fitting random forests for propensity",
  {
    fit <-
      chemical %>%
      ps_rf(poisox ~ age + sex, data = .)
    expect_s3_class(fit, "propmod")
    expect_s3_class(fit$model, "randomForest")
    expect_length(fit, 3)
  }
)
