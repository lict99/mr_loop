library(data.table)
library(openxlsx2)
library(future)

source("functions/mk_mr_loop_params.R", local = TRUE)
source("functions/run_mr_loop.R", local = TRUE)

plan(multisession)

options(future.globals.maxSize = 2000 * 1024^2)

output_dir <- "results/06"
dir.create(output_dir, FALSE, TRUE)

fmt_exp_args <- mk_fmt_params(
  type = "exposure",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  id_col = "id.exposure",
  chr_col = "chr.exposure",
  pos_col = "pos.exposure"
)

fmt_otc_args <- mk_fmt_params(
  type = "outcome",
  snp_col = "variant",
  beta_col = "coef",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "allele_eff",
  other_allele_col = "allele_ref",
  pval_col = "p_value",
  chr_col = "chr",
  pos_col = "pos"
)

clump_args <- mk_clump_params(
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EAS",
  bfile = "data/1kg.v3/EAS",
  plink_bin = Sys.getenv("PLINK_BIN")
)

results <- run_mr_loop(
  exposure_files = list.files("data/open_gwas/data", "csv$", full.names = TRUE),
  outcome_files = list.files("results/04", "results\\.csv$", full.names = TRUE),
  fmt_exposure_params = fmt_exp_args,
  fmt_outcome_params = fmt_otc_args,
  clump_params = clump_args
)

write_xlsx(results$mr_df, file.path(output_dir, "mr_loop_results.xlsx"))
write_xlsx(results$error_df, file.path(output_dir, "mr_loop_error.xlsx"))

fwrite(results$mr_df, file.path(output_dir, "mr_loop_results.csv"))
fwrite(results$error_df, file.path(output_dir, "mr_loop_error.csv"))
