# %% Env settings
library(qqman)
library(data.table)

output_dir <- "results/05"
dir.create(output_dir, FALSE, TRUE)

# %% Reading GWAS results and plotting
for (event in c("os", "css", "dfs")) {
  message(paste0("Reading ", event |> toupper(), " GWAS results"))
  gwas <- fread(
    grep(
      paste0("gwas_", event, ".+results\\.csv$"),
      list.files("results/04", full.names = TRUE),
      value = TRUE
    ),
    data.table = FALSE,
    na.strings = c("NA", ""),
    strip.white = TRUE,
    showProgress = FALSE
  )

  message("Plotting ", event |> toupper(), " GWAS results")
  png(
    file.path(output_dir, paste0("gwas_", event, "_manhattan.png")),
    res = 300,
    width = 8,
    height = 6,
    units = "in"
  )
  # Manhattan plot
  manhattan(
    gwas,
    chr = "chr",
    bp = "pos",
    p = "p_value",
    snp = "variant",
    annotatePval = 5e-8,
    annotateTop = FALSE,
    main = paste0(event |> toupper(), " GWAS Manhattan plot")
  )
  dev.off()
  message(paste0(event |> toupper(), " GWAS Manhattan plot saved"))

  message("Plotting ", event |> toupper(), " GWAS QQ plot")
  png(
    file.path(output_dir, paste0("gwas_", event, "_qq.png")),
    res = 300,
    width = 6,
    height = 6,
    units = "in"
  )
  # QQ plot
  qq(gwas[["p_value"]], main = paste0(event |> toupper(), " GWAS QQ plot"))
  dev.off()
  message(paste0(event |> toupper(), " GWAS QQ plot saved"))

  rm(gwas)
}
