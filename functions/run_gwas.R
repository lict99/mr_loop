run_gwas <- function(
    bed_file,
    model_params,
    output_file,
    chunk_size = 10000,
    min_maf = NULL) {
  loadNamespace("GWASTools")
  loadNamespace("SNPRelate")
  loadNamespace("future.apply")
  loadNamespace("data.table")
  loadNamespace("progress")

  if (!("model_params" %in% class(model_params))) {
    stop(
      "`model_params` must be an object of class 'model_params'",
      call. = FALSE
    )
  }

  gds_file <- tempfile(fileext = ".gds")

  SNPRelate::snpgdsBED2GDS(
    bed.fn = bed_file,
    out.gdsfn = gds_file
  )
  on.exit(unlink(gds_file), add = TRUE)

  gds <- GWASTools::GdsGenotypeReader(
    gds_file,
    YchromCode = 24L,
    XYchromCode = 25L
  )
  on.exit(GWASTools::close(gds), add = TRUE, after = FALSE)

  allele_a <- GWASTools::getAlleleA(gds)
  allele_b <- GWASTools::getAlleleB(gds)
  chromosome <- GWASTools::getChromosome(gds)
  position <- GWASTools::getPosition(gds)
  snp_id <- GWASTools::getSnpID(gds)
  n_snp <- GWASTools::nsnp(gds)

  scan_id <- GWASTools::getScanID(gds)

  covar_id <- model_params[["id"]]
  valid_id <- intersect(covar_id, scan_id)

  covar_id_idx <- match(valid_id, covar_id, 0)
  model_params[["id"]] <- covar_id[covar_id_idx]
  model_params[["data"]] <- model_params[["data"]][covar_id_idx, , drop = FALSE]
  gt_id_idx <- match(valid_id, scan_id, 0)
  if (!identical(model_params[["id"]], scan_id[gt_id_idx])) {
    stop("ID mismatch", call. = FALSE)
  }

  n_chunks <- ceiling(n_snp / chunk_size)

  pb <- progress::progress_bar$new(
    format = "Running :current/:total [:bar] :percent in :elapsedfull ETA :eta",
    total = n_chunks,
    clear = FALSE,
    show_after = 0
  )

  header <- c(
    "variant", "chr", "pos", "allele_eff", "allele_ref", "eaf", "maf"
  )

  invalid_file <- paste(output_file, "remove.csv", sep = "_")
  results_file <- paste(output_file, "results.csv", sep = "_")

  data.table::fwrite(
    as.data.frame(t(header)),
    invalid_file,
    col.names = FALSE
  )

  data.table::fwrite(
    as.data.frame(t(c(header, model_params[["header"]]))),
    results_file,
    col.names = FALSE
  )

  for (i in seq_len(n_chunks)) {
    pb$tick()

    if (i == n_chunks) {
      gt_start <- (i - 1) * chunk_size + 1
      gt_step <- -1
      gt_idx <- seq(gt_start, n_snp)
    } else {
      gt_start <- (i - 1) * chunk_size + 1
      gt_step <- chunk_size
      gt_idx <- seq(gt_start, gt_start + chunk_size - 1)
    }

    gt_a <- GWASTools::getGenotype(
      gds,
      snp = c(gt_start, gt_step),
      scan = c(1, -1),
      drop = FALSE,
      use.names = FALSE
    )
    gt_a <- gt_a[, gt_id_idx, drop = FALSE]

    fit_list <- future.apply::future_apply(
      gt_a, 1,
      function(gt) {
        gt_nm <- model_params[["gt_nm"]]
        covar_df <- model_params[["data"]]
        covar_df[[gt_nm]] <- gt
        fit_df <- na.omit(covar_df)
        eaf <- sum(fit_df[[gt_nm]]) / (2 * nrow(fit_df))
        maf <- if (eaf < 0.5) eaf else 1 - eaf
        is_valid_maf <- if (is.null(min_maf)) {
          if (maf >= 0 & maf <= 0.5) TRUE else FALSE
        } else {
          if (maf > min_maf & maf >= 0 & maf <= 0.5) TRUE else FALSE
        }
        if (!is_valid_maf) {
          return(data.frame(eaf = eaf, maf = maf))
        }

        fun <- model_params[["fun"]]
        fun_params <- model_params[["fun_params"]]
        fun_params[["data"]] <- fit_df
        model_df <- do.call(fun, fun_params)

        return(cbind(data.frame(eaf = eaf, maf = maf), model_df))
      },
      simplify = FALSE
    )

    is_valid_fit <- vapply(
      fit_list,
      function(x) if (ncol(x) == 2L) FALSE else TRUE,
      FUN.VALUE = logical(1L),
      USE.NAMES = FALSE
    )

    gt_info <- data.frame(
      variant = snp_id[gt_idx],
      chr = chromosome[gt_idx],
      pos = position[gt_idx],
      allele_eff = allele_a[gt_idx],
      allele_ref = allele_b[gt_idx]
    )

    invalid_gt <- cbind(
      gt_info[!is_valid_fit, , drop = FALSE],
      do.call(rbind, fit_list[!is_valid_fit])
    )

    model_results <- cbind(
      gt_info[is_valid_fit, , drop = FALSE],
      do.call(rbind, fit_list[is_valid_fit])
    )

    data.table::fwrite(invalid_gt, invalid_file, append = TRUE)
    data.table::fwrite(model_results, results_file, append = TRUE)
  }
  return(invisible(NULL))
}
