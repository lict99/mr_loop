.fmt_exp_data <- function(
    exposure_files,
    fmt_exposure_params) {
  message("Reading exposure data...")
  fmt_exp_data <- lapply(
    exposure_files,
    function(file) {
      data <- data.table::fread(
        file,
        data.table = FALSE,
        showProgress = FALSE
      )

      id_col_arg <- fmt_exposure_params[["id_col"]]
      if (is.null(id_col_arg) || !(id_col_arg %in% colnames(data))) {
        id_col_dft <- as.list(formals(TwoSampleMR::format_data))$id_col
        fmt_exposure_params[["id_col"]] <- id_col_dft
        if (is.null(id_col_dft) || !(id_col_dft %in% colnames(data))) {
          data[["id"]] <- basename(file)
          fmt_exposure_params[["id_col"]] <- "id"
        }
      }

      phe_col_arg <- fmt_exposure_params[["phenotype_col"]]
      if (is.null(phe_col_arg) || !(phe_col_arg %in% colnames(data))) {
        phe_col_dft <- as.list(formals(TwoSampleMR::format_data))$phenotype_col
        fmt_exposure_params[["phenotype_col"]] <- phe_col_dft
        if (is.null(phe_col_dft) || !(phe_col_dft %in% colnames(data))) {
          data[["phenotype"]] <- basename(file)
          fmt_exposure_params[["phenotype_col"]] <- "phenotype"
        }
      }

      fmt_exposure_params[["dat"]] <- data

      fmt_data <- tryCatch(
        do.call(TwoSampleMR::format_data, fmt_exposure_params),
        error = function(e) {
          id_exp <- data[[fmt_exposure_params[["id_col"]]]][[1]]
          return(
            data.frame(
              id = if (is.na(id_exp)) basename(file) else id_exp,
              info = paste0(
                "When formatting exposure data with file `",
                basename(file), "`, ",
                trimws(as.character(e))
              )
            )
          )
        }
      )

      return(fmt_data)
    }
  )

  is_valid_exp <- vapply(
    fmt_exp_data,
    function(x) !(ncol(x) == 2L && "info" %in% colnames(x)),
    FUN.VALUE = logical(1L),
    USE.NAMES = FALSE
  )

  error_df_exp <- do.call(rbind, fmt_exp_data[!is_valid_exp])
  fmt_exp_data <- fmt_exp_data[is_valid_exp]

  return(list(data = fmt_exp_data, error_df = error_df_exp))
}

.fmt_otc_data <- function(
    outcome_files,
    fmt_outcome_params) {
  message("Reading outcome data...")
  fmt_otc_data <- lapply(
    outcome_files,
    function(file) {
      data <- data.table::fread(
        file,
        data.table = FALSE,
        showProgress = FALSE
      )

      id_col_arg <- fmt_outcome_params[["id_col"]]
      if (is.null(id_col_arg) || !(id_col_arg %in% colnames(data))) {
        id_col_dft <- as.list(formals(TwoSampleMR::format_data))$id_col
        fmt_outcome_params[["id_col"]] <- id_col_dft
        if (is.null(id_col_dft) || !(id_col_dft %in% colnames(data))) {
          data[["id"]] <- basename(file)
          fmt_outcome_params[["id_col"]] <- "id"
        }
      }

      phe_col_arg <- fmt_outcome_params[["phenotype_col"]]
      if (is.null(phe_col_arg) || !(phe_col_arg %in% colnames(data))) {
        phe_col_dft <- as.list(formals(TwoSampleMR::format_data))$phenotype_col
        fmt_outcome_params[["phenotype_col"]] <- phe_col_dft
        if (is.null(phe_col_dft) || !(phe_col_dft %in% colnames(data))) {
          data[["phenotype"]] <- basename(file)
          fmt_outcome_params[["phenotype_col"]] <- "phenotype"
        }
      }

      fmt_outcome_params[["dat"]] <- data

      fmt_data <- tryCatch(
        do.call(TwoSampleMR::format_data, fmt_outcome_params),
        error = function(e) {
          id_otc <- data[[fmt_outcome_params[["id_col"]]]][[1]]
          return(
            data.frame(
              id = if (is.na(id_otc)) basename(file) else id_otc,
              info = paste0(
                "When formatting outcome data with file `",
                basename(file), "`, ",
                trimws(as.character(e))
              )
            )
          )
        }
      )

      return(fmt_data)
    }
  )

  is_valid_otc <- vapply(
    fmt_otc_data,
    function(x) !(ncol(x) == 2L && "info" %in% colnames(x)),
    FUN.VALUE = logical(1L),
    USE.NAMES = FALSE
  )

  error_df_otc <- do.call(rbind, fmt_otc_data[!is_valid_otc])
  fmt_otc_data <- fmt_otc_data[is_valid_otc]

  return(list(data = fmt_otc_data, error_df = error_df_otc))
}

.clump_exp_data <- function(
    fmt_exp_data,
    exposure_p,
    clump_params) {
  message("Clumping exposure data...")

  clumped_exp_data <- lapply(
    fmt_exp_data,
    function(df) {
      if (!("pval.exposure" %in% colnames(df))) {
        stop("No p-value column in the exposure data", call. = FALSE)
      }
      if (!is.numeric(df[["pval.exposure"]])) {
        stop(
          "P-value column is not numeric in the exposure data",
          call. = FALSE
        )
      }

      df_filter <- df[df[["pval.exposure"]] <= exposure_p, , drop = FALSE]

      if (nrow(df_filter) == 0L) {
        return(
          data.frame(
            id = df[["id.exposure"]][[1]],
            info = "No significant exposure SNP with below given p-value"
          )
        )
      }

      clump_params[["dat"]] <- df_filter
      clump_df <- tryCatch(
        do.call(TwoSampleMR::clump_data, clump_params),
        error = function(e) {
          data.frame(
            id = df[["id.exposure"]][[1]],
            info = paste(
              "When clumping exposure data,",
              trimws(as.character(e))
            )
          )
        }
      )

      if (nrow(clump_df) == 0L) {
        return(
          data.frame(
            id = df[["id.exposure"]][[1]],
            info = "No exposure SNP after clumping"
          )
        )
      }

      return(clump_df)
    }
  )

  is_valid_clump <- vapply(
    clumped_exp_data,
    function(x) !(ncol(x) == 2L && ("info" %in% colnames(x))),
    FUN.VALUE = logical(1L),
    USE.NAMES = FALSE
  )

  error_df_clump <- do.call(rbind, clumped_exp_data[!is_valid_clump])
  clumped_exp_data <- clumped_exp_data[is_valid_clump]

  return(list(data = clumped_exp_data, error_df = error_df_clump))
}

.harmonise_data <- function(
    clumped_exp_data,
    fmt_otc_data,
    harmonise_params) {
  message("Harmonising data...")

  combined_idx <- expand.grid(
    exposure = seq_along(clumped_exp_data),
    outcome = seq_along(fmt_otc_data)
  )

  harmonise_data <- lapply(
    seq_len(nrow(combined_idx)),
    function(idx) {
      exp_idx <- combined_idx[idx, ][[1]]
      otc_idx <- combined_idx[idx, ][[2]]
      exp_df <- clumped_exp_data[[exp_idx]]
      otc_df <- fmt_otc_data[[otc_idx]]
      if (length(intersect(exp_df[["SNP"]], otc_df[["SNP"]])) == 0L) {
        return(
          data.frame(
            id = paste(
              exp_df[["id.exposure"]][[1]],
              otc_df[["id.outcome"]][[1]],
              sep = " | "
            ),
            info = "No common SNP when harmonising"
          )
        )
      }
      harmonise_params[["exposure_dat"]] <- exp_df
      harmonise_params[["outcome_dat"]] <- otc_df
      harmonise_df <- tryCatch(
        do.call(TwoSampleMR::harmonise_data, harmonise_params),
        error = function(e) {
          data.frame(
            id = paste(
              exp_df[["id.exposure"]][[1]],
              otc_df[["id.outcome"]][[1]],
              sep = " | "
            ),
            info = paste("When harmonising data,", trimws(as.character(e)))
          )
        }
      )

      if (sum(harmonise_df[["mr_keep"]]) == 0L) {
        harmonise_df <- data.frame(
          id = paste(
            exp_df[["id.exposure"]][[1]],
            otc_df[["id.outcome"]][[1]],
            sep = " | "
          ),
          info = "No SNP available for MR analysis after harmonising"
        )
      }

      return(harmonise_df)
    }
  )

  is_valid_harmonise <- vapply(
    harmonise_data,
    function(x) !(ncol(x) == 2L && "info" %in% colnames(x)),
    FUN.VALUE = logical(1L),
    USE.NAMES = FALSE
  )

  error_df_harmonise <- do.call(rbind, harmonise_data[!is_valid_harmonise])
  harmonise_data <- harmonise_data[is_valid_harmonise]

  return(list(data = harmonise_data, error_df = error_df_harmonise))
}

.run_mr <- function(
    harmonise_data,
    mr_params) {
  message("Running MR analysis...")
  mr_data <- future.apply::future_lapply(
    harmonise_data,
    function(df) {
      tryCatch(
        {
          set.seed(1)
          mr_params[["dat"]] <- df
          mr_res <- do.call(TwoSampleMR::mr, mr_params)
          TwoSampleMR::generate_odds_ratios(mr_res)
        },
        error = function(e) {
          data.frame(
            id = paste(
              df[["id.exposure"]][[1]],
              df[["id.outcome"]][[1]],
              sep = " | "
            ),
            info = paste("When Running MR,", trimws(as.character(e)))
          )
        }
      )
    },
    future.packages = c("TwoSampleMR"),
    future.seed = TRUE
  )

  is_valid_mr <- vapply(
    mr_data,
    function(x) !(ncol(x) == 2L && "info" %in% colnames(x)),
    FUN.VALUE = logical(1L),
    USE.NAMES = FALSE
  )

  error_df_mr <- do.call(rbind, mr_data[!is_valid_mr])
  mr_data <- mr_data[is_valid_mr]

  return(list(data = mr_data, error_df = error_df_mr))
}

run_mr_loop <- function(
    exposure_files,
    outcome_files,
    exposure_p = 5e-8,
    fmt_exposure_params = mk_fmt_params(), # nolint: object_usage_linter.
    fmt_outcome_params = mk_fmt_params(), # nolint: object_usage_linter.
    clump_params = mk_clump_params(), # nolint: object_usage_linter.
    harmonise_params = mk_harmonise_params(), # nolint: object_usage_linter.
    mr_params = mk_mr_params()) { # nolint: object_usage_linter.
  loadNamespace("TwoSampleMR")
  loadNamespace("data.table")
  loadNamespace("future.apply")

  fmt_exp_data <- .fmt_exp_data(
    exposure_files = exposure_files,
    fmt_exposure_params = fmt_exposure_params
  )

  error_df_exp <- fmt_exp_data[["error_df"]]
  fmt_exp_data <- fmt_exp_data[["data"]]

  fmt_otc_data <- .fmt_otc_data(
    outcome_files = outcome_files,
    fmt_outcome_params = fmt_outcome_params
  )

  error_df_otc <- fmt_otc_data[["error_df"]]
  fmt_otc_data <- fmt_otc_data[["data"]]

  clumped_exp_data <- .clump_exp_data(
    fmt_exp_data = fmt_exp_data,
    exposure_p = exposure_p,
    clump_params = clump_params
  )

  error_df_clump <- clumped_exp_data[["error_df"]]
  clumped_exp_data <- clumped_exp_data[["data"]]

  harmonise_data <- .harmonise_data(
    clumped_exp_data = clumped_exp_data,
    fmt_otc_data = fmt_otc_data,
    harmonise_params = harmonise_params
  )

  error_df_harmonise <- harmonise_data[["error_df"]]
  harmonise_data <- harmonise_data[["data"]]

  mr_data <- .run_mr(
    harmonise_data = harmonise_data,
    mr_params = mr_params
  )

  error_df_mr <- mr_data[["error_df"]]
  mr_data <- mr_data[["data"]]

  mr_df <- do.call(rbind, mr_data)

  error_df <- rbind(
    error_df_otc,
    error_df_exp,
    error_df_clump,
    error_df_harmonise,
    error_df_mr
  )

  return(list(mr_df = mr_df, error_df = error_df))
}
