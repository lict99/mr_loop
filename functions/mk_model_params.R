.chk_model_params <- function(ignore, ...) {
  input_params <- list(...)
  if (length(input_params) > 0 && is.null(names(input_params))) {
    stop("Do not support positional arguments in `...`", call. = FALSE)
  }
  if (length(input_params) > 0 && any("" %in% names(input_params))) {
    stop("Do not support positional arguments in `...`", call. = FALSE)
  }
  for (i in ignore) {
    if (i %in% names(input_params)) {
      stop("Parameter `", i, "` should be ignored", call. = FALSE)
    }
  }
  return(TRUE)
}

mk_coxph_params <- function(
    covariates_df,
    id_col,
    time_col,
    event_col,
    covariates_col,
    ...) {
  loadNamespace("survival")
  extract_coxph <- function(formula, data, ...) {
    fit <- survival::coxph(formula = formula, data = data, ...)
    fit_smr <- summary(fit)
    return(
      data.frame(
        coef = fit_smr[["coefficients"]]["genotype", "coef"],
        se = fit_smr[["coefficients"]]["genotype", "se(coef)"],
        p_value = fit_smr[["coefficients"]]["genotype", "Pr(>|z|)"],
        hr = fit_smr[["conf.int"]]["genotype", "exp(coef)"],
        hr_l95 = fit_smr[["conf.int"]]["genotype", "lower .95"],
        hr_u95 = fit_smr[["conf.int"]]["genotype", "upper .95"],
        n = fit_smr[["n"]],
        n_event = fit_smr[["nevent"]]
      )
    )
  }

  gt_nm <- "genotype"

  .chk_model_params(c("formula", "data"), ...)
  dot_params <- list(...)
  data <- covariates_df[c(time_col, event_col, covariates_col)]
  y_char <- sprintf(
    "survival::Surv(time = %s, event = %s)",
    time_col, event_col
  )
  x_char <- paste(c(gt_nm, covariates_col), collapse = " + ")
  fml <- as.formula(paste(y_char, x_char, sep = " ~ "))

  return(
    structure(
      list(
        id = covariates_df[[id_col]],
        gt_nm = gt_nm,
        data = data,
        fun = extract_coxph,
        fun_params = c(list(formula = fml, data = NULL), dot_params),
        header = c(
          "coef", "se", "p_value", "hr", "hr_l95", "hr_u95", "n", "n_event"
        )
      ),
      class = "model_params"
    )
  )
}
