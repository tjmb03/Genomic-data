# R/reporting_tables.R
make_benchmark_table <- function(metrics_summary, runtimes = NULL) {
  df <- metrics_summary
  if (!is.null(runtimes)) {
    df <- dplyr::left_join(df, runtimes, by = "method")
  } else {
    df$runtime_sec <- NA_real_
  }

  df <- df %>%
    mutate(
      batch_mixing = case_when(
        is.na(batch_sil_mean) ~ NA_character_,
        batch_sil_mean <= 0.05 ~ "Excellent",
        batch_sil_mean <= 0.15 ~ "Good",
        batch_sil_mean <= 0.30 ~ "Moderate",
        TRUE ~ "Poor"
      ),
      batch_effect_r2 = case_when(
        is.na(batch_r2_mean) ~ NA_character_,
        batch_r2_mean <= 0.05 ~ "Low",
        batch_r2_mean <= 0.15 ~ "Moderate",
        TRUE ~ "High"
      ),
      biology_conservation = case_when(
        is.na(cluster_sil_mean) ~ NA_character_,
        cluster_sil_mean >= 0.35 ~ "Excellent",
        cluster_sil_mean >= 0.20 ~ "Good",
        cluster_sil_mean >= 0.10 ~ "Moderate",
        TRUE ~ "Weak"
      ),
      biology_signal_r2 = case_when(
        is.na(cluster_r2_mean) ~ NA_character_,
        cluster_r2_mean >= 0.20 ~ "High",
        cluster_r2_mean >= 0.10 ~ "Moderate",
        TRUE ~ "Low"
      )
    ) %>%
    select(method, reduction_used, n_cells_total, n_cells_metrics, dims_used,
           batch_sil_mean, batch_mixing,
           batch_r2_mean, batch_effect_r2,
           cluster_sil_mean, biology_conservation,
           cluster_r2_mean, biology_signal_r2,
           runtime_sec)

  df
}
