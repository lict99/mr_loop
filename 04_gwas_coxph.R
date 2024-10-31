library(future)
library(fastDummies)

if (supportsMulticore()) plan(multicore) else plan(multisession)

options(future.globals.maxSize = 1000 * 1024^2)

source("functions/mk_model_params.R", local = TRUE)
source("functions/run_gwas.R", local = TRUE)

output_dir <- "results/04"

dir.create(output_dir, FALSE, TRUE)

load("results/01/hx_df.rda")

pca <- read.delim(
  "results/03/pca5.eigenvec",
  header = FALSE,
  sep = " "
) |>
  setNames(c("FID", "IID", paste0("PC", 1:5)))

covariates_df <- merge(
  pca,
  hx_df,
  by.x = "IID",
  by.y = "id",
  all.x = TRUE,
  sort = FALSE
) |>
  dummy_cols(
    select_columns = "sex",
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE,
    ignore_na = TRUE
  )

for (event in c("os", "css", "dfs")) {
  time <- paste0(event, "_time_days")
  covariates <- c("age", grep("PC|sex", names(covariates_df), value = TRUE))
  df <- covariates_df[c("IID", covariates, time, event)]
  df <- na.omit(df)
  df <- df[df[[time]] > 30, , drop = FALSE]

  coxph_params <- mk_coxph_params(
    covariates_df = df,
    id_col = "IID",
    time_col = time,
    event_col = event,
    covariates_col = covariates,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )

  run_gwas(
    bed_file = "results/03/data_final.bed",
    model_params = coxph_params,
    output_file = file.path(
      output_dir,
      paste("gwas", event, "coxph", sep = "_")
    ),
    chunk_size = 10000,
    min_maf = NULL
  )
}
