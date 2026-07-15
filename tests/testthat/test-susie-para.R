test_that("SuSiE defaults preserve package settings and native omissions", {
  para <- SuSiEIRLS:::.susie_default_para()

  expect_identical(
    para,
    list(
      scaled_prior_variance = 2,
      estimate_residual_variance = TRUE,
      residual_variance = 0.5,
      residual_variance_lowerbound = 0.1,
      residual_variance_upperbound = 1.01,
      estimate_prior_variance = TRUE,
      estimate_prior_method = "optim",
      max_iter = 30,
      coverage = 0.9
    )
  )
  expect_false("min_abs_corr" %in% names(para))
  expect_false("check_null_threshold" %in% names(para))
})

test_that("susie_para supports partial native overrides", {
  user_para <- SuSiEIRLS:::.resolve_susie_para(list(
    estimate_residual_variance = FALSE,
    min_abs_corr = 0.2,
    max_iter = 7
  ))
  structural <- list(
    XtX = diag(2), Xty = c(1, 0), yty = 1, n = 10, L = 2
  )
  para <- SuSiEIRLS:::.susie_iteration_args(
    user_para, structural, iter = 3, min.iter = 2
  )

  expect_false(para$estimate_residual_variance)
  expect_equal(para$min_abs_corr, 0.2)
  expect_equal(para$max_iter, 7)
  expect_equal(para$coverage, 0.9)
})

test_that("susie_para rejects structural and invalid arguments", {
  expect_error(
    SuSiEIRLS:::.resolve_susie_para(list(XtX = diag(2))),
    "cannot set structural"
  )
  expect_error(
    SuSiEIRLS:::.resolve_susie_para(list(L = 2)),
    "cannot set structural"
  )
  expect_error(
    SuSiEIRLS:::.resolve_susie_para(list(not_a_susie_argument = 1)),
    "Unknown susieR::susie_ss parameter"
  )
  expect_error(
    SuSiEIRLS:::.resolve_susie_para(list(1)),
    "named list"
  )
})

test_that("iteration controls retain warm-up behavior", {
  para <- SuSiEIRLS:::.resolve_susie_para(list(
    scaled_prior_variance = 1.5,
    estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE,
    max_iter = 7,
    coverage = 0.8
  ))
  structural <- list(
    XtX = diag(2), Xty = c(1, 0), yty = 1, n = 10, L = 2
  )

  warm <- SuSiEIRLS:::.susie_iteration_args(para, structural, 2, 2)
  update <- SuSiEIRLS:::.susie_iteration_args(para, structural, 3, 2)
  expect_equal(warm$scaled_prior_variance, 2)
  expect_false(warm$estimate_prior_variance)
  expect_true(warm$estimate_residual_variance)
  expect_equal(warm$max_iter, 30)
  expect_equal(warm$coverage, 0.9)
  expect_equal(update$scaled_prior_variance, 1.5)
  expect_false(update$estimate_prior_variance)
  expect_false(update$estimate_residual_variance)
  expect_equal(update$max_iter, 7)
  expect_equal(update$coverage, 0.8)
  expect_identical(update$XtX, structural$XtX)

  default_update <- SuSiEIRLS:::.susie_iteration_args(
    SuSiEIRLS:::.resolve_susie_para(), structural, 3, 2
  )
  expect_equal(default_update$scaled_prior_variance, 2)
  expect_true(default_update$estimate_prior_variance)
})

test_that("outer iteration default is 10", {
  expect_equal(formals(SuSiE_IRLS)$max.iter, 10)
})

test_that("susie_para is the only SuSiE parameter interface", {
  expect_false("..." %in% names(formals(SuSiE_IRLS)))
  expect_error(
    SuSiEIRLS:::.resolve_susie_para(list(susie.iter = 20)),
    "Unknown susieR::susie_ss parameter"
  )
})
