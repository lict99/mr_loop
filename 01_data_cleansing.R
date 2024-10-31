# %% Env settings
source("functions/char2date.R", local = TRUE)

output_dir <- "results/01"
dir.create(output_dir, FALSE, TRUE)

# %% Reading raw data with comments
raw_with_cmmts <- read.csv(
  file = "data/raw_data_with_comments_240918.csv",
  na.strings = c("", " ", "NA", "#N/A"),
  strip.white = TRUE
)

# %% Arrange data to correct format
hx_df <- with(
  # We only include valid patients.
  raw_with_cmmts[raw_with_cmmts$inclusion == 1L, , drop = FALSE],
  data.frame(
    name = as.character(name),
    id = as.character(id),
    register = as.character(register),
    sex = as.integer(sex) |>
      factor(c(2L, 1L), c("female", "male")) |>
      as.character(),
    age = as.numeric(age),
    surgery_date = char2date(surgery_date),
    diagnosis_date_for_nonsurgery = char2date(diagnosis_date_for_nonsurgery),
    c_stage = as.character(c_stage),
    p_stage = as.character(p_stage),
    death = as.integer(death),
    death_by_crc = as.integer(death_by_crc),
    last_fu_date = char2date(last_fu_date),
    death_date = char2date(death_date),
    local_recurrence = as.integer(local_recurrence),
    local_recurrence_date = char2date(local_recurrence_date),
    metastasis = as.integer(metastasis),
    metastasis_date = char2date(metastasis_date),
    last_radio_date = char2date(last_radio_date)
  )
)

# %% OS indicator and time
hx_df[["os"]] <- with(hx_df, death)
hx_df[["os_time_days"]] <- with(
  hx_df,
  {
    # OS start date is primarily surgery date if available, otherwise diagnosis
    # date.
    os_start_date <- ifelse(
      !is.na(surgery_date),
      surgery_date |> as.numeric(),
      diagnosis_date_for_nonsurgery |> as.numeric()
    )
    # OS end date is death date for deceased patients, otherwise the latest date
    # of follow-up or radiography test for living patients.
    os_end_date <- ifelse(
      os == 1L,
      death_date |> as.numeric(),
      ifelse(
        os == 0L,
        pmax(
          last_fu_date |> as.numeric(),
          last_radio_date |> as.numeric(),
          na.rm = TRUE
        ),
        NA_real_
      )
    )
    # OS time is the difference between OS end date and OS start date.
    os_end_date - os_start_date
  }
)

# %% CSS indicator and time
hx_df[["css"]] <- with(hx_df, death_by_crc)
hx_df[["css_time_days"]] <- with(
  hx_df,
  {
    # CSS start date is primarily surgery date if available, otherwise diagnosis
    # date.
    css_start_date <- ifelse(
      !is.na(surgery_date),
      surgery_date |> as.numeric(),
      diagnosis_date_for_nonsurgery |> as.numeric()
    )
    # CSS end date is death date for patients who died of CRC, otherwise the
    # latest date of death, follow-up or radiography test for patients without
    # CSS events.
    css_end_date <- ifelse(
      css == 1L,
      death_date |> as.numeric(),
      ifelse(
        css == 0L,
        pmax(
          last_fu_date |> as.numeric(),
          last_radio_date |> as.numeric(),
          death_date |> as.numeric(),
          na.rm = TRUE
        ),
        NA_real_
      )
    )
    # CSS time is the difference between CSS end date and CSS start date.
    css_end_date - css_start_date
  }
)

# %% DFS indicator and time
hx_df[["dfs"]] <- with(
  hx_df,
  ifelse(
    # We only consider patients who had surgery.
    !is.na(surgery_date),
    ifelse(
      death_by_crc == 1L | local_recurrence == 1L | metastasis == 1L,
      1L,
      ifelse(
        death_by_crc == 0L & local_recurrence == 0L & metastasis == 0L,
        0L,
        NA_integer_
      )
    ),
    NA_integer_
  )
)
hx_df[["dfs_time_days"]] <- with(
  hx_df,
  {
    # DFS start date is surgery date
    dfs_start_date <- surgery_date |> as.numeric()
    # DFS end date is earliest date of death by CRC, local recurrence, or
    # metastasis for patients who had events, otherwise the latest date of
    # radiography test for patients without events.
    # We do not consider the last follow-up date for patients without DFS events
    # because we cannot detect DFS events according to the phone contact
    # follow-up.
    dfs_end_date <- ifelse(
      dfs == 1L,
      pmin(
        death_date |> as.numeric(),
        local_recurrence_date |> as.numeric(),
        metastasis_date |> as.numeric(),
        na.rm = TRUE
      ),
      ifelse(
        dfs == 0L,
        last_radio_date |> as.numeric(),
        NA_real_
      )
    )
    # DFS time is the difference between DFS end date and DFS start date.
    dfs_end_date - dfs_start_date
  }
)

# %% Saving data
save(hx_df, file = file.path(output_dir, "hx_df.rda"))
