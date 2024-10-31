# %% Env settings
load("results/01/hx_df.rda")

output_dir <- "results/02"
dir.create(output_dir, FALSE, TRUE)

# %% Reading samples of SNP array
snp_samples <- read.delim(
  "data/SNParray/all/raw_data/22B0218D_23B0203E_23P1212B00A_24B0516F.fam",
  sep = " ",
  header = FALSE
)

# %% Getting shared samples
shared_samples <- intersect(
  snp_samples$V1,
  hx_df$id
)

# %% Getting sex information for shared samples
sex_info <- hx_df[match(shared_samples, hx_df$id), c("id", "id", "sex")] |>
  transform(
    sex = factor(sex, levels = c("male", "female"), labels = c("1", "2")) |>
      as.character()
  )

# %% Saving sex information for PLINK
write.table(
  sex_info,
  file.path(output_dir, "sex_info.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  na = "0"
)
