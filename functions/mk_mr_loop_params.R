.mk_params <- function(fun, ignore, ...) {
  default_params <- as.list(formals(match.fun(fun)))

  if ("..." %in% names(default_params)) {
    stop("Do not support function with `...` parameter", call. = FALSE)
  }

  input_params <- list(...)

  for (i in ignore) {
    if (i %in% names(input_params)) {
      stop("Parameter `", i, "` should be ignored", call. = FALSE)
    }
  }

  if (length(input_params) == 0) {
    return(input_params)
  }

  if (length(input_params) > 0 && is.null(names(input_params))) {
    stop("Do not support positional arguments", call. = FALSE)
  }

  if (length(input_params) > 0 && any("" %in% names(input_params))) {
    stop("Do not support positional arguments", call. = FALSE)
  }

  if (!all(names(input_params) %in% names(default_params))) {
    stop("Do not support unknown parameters", call. = FALSE)
  }

  valid_input_param_nms <- intersect(names(input_params), names(default_params))
  params <- input_params[valid_input_param_nms]

  return(params)
}

mk_fmt_params <- function(...) {
  params <- .mk_params(TwoSampleMR::format_data, c("dat"), ...)
  class(params) <- "fmt_params"
  return(params)
}

mk_clump_params <- function(...) {
  params <- .mk_params(TwoSampleMR::clump_data, c("dat"), ...)
  class(params) <- "clump_params"
  return(params)
}

mk_harmonise_params <- function(...) {
  params <- .mk_params(
    TwoSampleMR::harmonise_data,
    c("exposure_dat", "outcome_dat"),
    ...
  )
  class(params) <- "harmonise_params"
  return(params)
}

mk_mr_params <- function(...) {
  params <- .mk_params(TwoSampleMR::mr, c("dat"), ...)
  class(params) <- "mr_params"
  return(params)
}
