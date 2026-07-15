test_that("susie_para preserves package defaults and native omissions", {
  para <- SuSiEIRLS:::.resolve_susie_para()

  expect_identical(
    para,
    list(
      estimate_residual_variance = TRUE,
      residual_variance = 0.5,
      residual_variance_lowerbound = 0.1,
      residual_variance_upperbound = 1,
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
  para <- SuSiEIRLS:::.resolve_susie_para(list(
    estimate_residual_variance = FALSE,
    min_abs_corr = 0.2,
    max_iter = 7
  ))

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
  para <- SuSiEIRLS:::.resolve_susie_para()
  structural <- list(
    XtX = diag(2), Xty = c(1, 0), yty = 1, n = 10, L = 2
  )

  warm <- SuSiEIRLS:::.susie_iteration_args(para, structural, 2, 2)
  update <- SuSiEIRLS:::.susie_iteration_args(para, structural, 3, 2)
  expect_equal(warm$scaled_prior_variance, 2)
  expect_false(warm$estimate_prior_variance)
  expect_equal(update$scaled_prior_variance, 3)
  expect_true(update$estimate_prior_variance)
  expect_identical(update$XtX, structural$XtX)

  fixed <- SuSiEIRLS:::.resolve_susie_para(list(
    scaled_prior_variance = 1.5,
    estimate_prior_variance = FALSE
  ))
  fixed <- SuSiEIRLS:::.susie_iteration_args(fixed, structural, 3, 2)
  expect_equal(fixed$scaled_prior_variance, 1.5)
  expect_false(fixed$estimate_prior_variance)
})

test_that("susie_para is the only SuSiE parameter interface", {
  expect_false("..." %in% names(formals(SuSiE_IRLS)))
  expect_error(
    SuSiEIRLS:::.resolve_susie_para(list(susie.iter = 20)),
    "Unknown susieR::susie_ss parameter"
  )
})
