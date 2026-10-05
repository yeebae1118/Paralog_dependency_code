# ATAC-seq TSS stabilized log2 fold-change analysis
# Input columns: group, bin, mean, sd, n, se, condition, gene_label, chr

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(mgcv)
library(ggplot2)
library(readr)

input_file <- "ATACseq_TSS_df.xlsx"
output_dir <- "ATACseq_TSS_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_excel(input_file, sheet = "Sheet1") %>%
  mutate(
    bin_num = as.numeric(bin),
    condition = factor(condition, levels = c("neutral", "loss")),
    gene_label = factor(gene_label, levels = c("nonsig", "sig"))
  )

# One global pseudocount is used for every chromosome arm and gene group so the
# comparisons remain on the same scale. The fifth percentile of all positive
# signals stabilizes ratios when either group has a very small mean signal.
positive_signals <- df$mean[is.finite(df$mean) & df$mean > 0]

if (length(positive_signals) == 0) {
  stop("No positive mean signals are available to estimate a pseudocount.")
}

pseudocount <- as.numeric(quantile(
  positive_signals,
  probs = 0.05,
  na.rm = TRUE,
  names = FALSE
))

write_csv(
  tibble(
    pseudocount_method = "5th percentile of all positive mean signals",
    pseudocount = pseudocount
  ),
  file.path(output_dir, "log2FC_pseudocount.csv")
)

# Confirm that all four groups exist for every chromosome and bin.
data_check <- df %>%
  count(chr, gene_label, condition, name = "rows") %>%
  complete(chr, gene_label, condition, fill = list(rows = 0))

write_csv(data_check, file.path(output_dir, "data_check.csv"))

# 1. Stabilized loss-versus-neutral log2 fold change for each paralog group.
# Negative values mean reduced accessibility in chromosome-loss samples;
# positive values mean increased accessibility in chromosome-loss samples.
difference_df <- df %>%
  select(chr, gene_label, bin_num, condition, mean, se) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(mean, se)
  ) %>%
  mutate(
    log2FC = log2(
      (mean_loss + pseudocount) /
        (mean_neutral + pseudocount)
    ),

    # Delta-method variance for log2((loss+p)/(neutral+p)), assuming the
    # summarized loss and neutral means are independent.
    variance_log2FC =
      (se_loss / ((mean_loss + pseudocount) * log(2)))^2 +
      (se_neutral / ((mean_neutral + pseudocount) * log(2)))^2,

    se_log2FC = sqrt(variance_log2FC),
    lower_log2FC_CI = log2FC - 1.96 * se_log2FC,
    upper_log2FC_CI = log2FC + 1.96 * se_log2FC,

    # Optional percentage-scale interpretation derived from the stabilized
    # ratio. Statistical tests below use log2FC, not percentage change.
    stabilized_percent_change = 100 * (2^log2FC - 1),
    low_signal = pmax(mean_loss, mean_neutral) < pseudocount,

    # Generic aliases retained so the downstream GAM/permutation code works.
    difference = log2FC,
    variance_difference = variance_log2FC,
    se_difference = sqrt(variance_difference),
    lower_difference = difference - 1.96 * se_difference,
    upper_difference = difference + 1.96 * se_difference
  )

write_csv(
  difference_df,
  file.path(output_dir, "stabilized_log2FC_curves.csv")
)


##### perform the pertumation test around TSS region 
# the point is we want to compare between sig and non-sig paralog genes
difference_df %>% 
  dplyr::filter(bin_num %in% c(55,56,57,58,59,60,
                               61,62,63,64,65)) %>% 
  dplyr::select(chr, gene_label, bin_num, log2FC)-> TSS_region_FC

# Perform permutation test 

paired_permutation_test <- function(group_A, group_B,
                                    n_perm = 10000,
                                    seed = 123) {
  stopifnot(length(group_A) == length(group_B))
  
  keep <- is.finite(group_A) & is.finite(group_B)
  difference <- group_A[keep] - group_B[keep]
  
  set.seed(seed)
  
  observed_difference <- mean(difference)
  
  permuted_differences <- replicate(n_perm, {
    signs <- sample(c(-1, 1), length(difference), replace = TRUE)
    mean(difference * signs)
  })
  
  p_value <- (
    sum(permuted_differences >= observed_difference) + 1
  ) / (n_perm + 1)
  
  list(
    observed_difference = observed_difference,
    p_value = p_value,
    n_pairs = length(difference),
    alternative = "Group A > Group B"
  )
}

result <- paired_permutation_test(
  group_A = TSS_region_FC[TSS_region_FC$gene_label == "nonsig",4]$log2FC,
  group_B = TSS_region_FC[TSS_region_FC$gene_label == "sig",4]$log2FC,
  n_perm = 100000
)


# -----------------------------------------------------------------------------
# OVERALL ANALYSIS: pool chromosome arms before inspecting individual arms
# -----------------------------------------------------------------------------
# Each chromosome arm is treated as one replicate. This prevents chromosome arms
# with more loss samples from dominating the overall profile.

upstream <- 3000
downstream <- 3000
total_bins <- max(df$bin_num, na.rm = TRUE)
up_bins <- total_bins * upstream / (upstream + downstream)
tss_pos <- up_bins + 1

# Overall stabilized log2FC curves. Each chromosome arm contributes one
# log2FC curve, so absolute signal differences among arms or gene
# groups do not dominate the combined comparison.
overall_log2FC <- difference_df %>%
  group_by(gene_label, bin_num) %>%
  summarise(
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    median_log2FC = median(log2FC, na.rm = TRUE),
    mean_stabilized_percent_change = mean(
      stabilized_percent_change,
      na.rm = TRUE
    ),
    n_arms = sum(!is.na(log2FC)),
    variance_across_arms = var(log2FC, na.rm = TRUE),
    se_across_arms = sqrt(variance_across_arms / n_arms),
    lower_CI = mean_log2FC -
      qt(0.975, df = pmax(n_arms - 1, 1)) * se_across_arms,
    upper_CI = mean_log2FC +
      qt(0.975, df = pmax(n_arms - 1, 1)) * se_across_arms,
    low_signal_arms = sum(low_signal, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  overall_log2FC,
  file.path(output_dir, "overall_stabilized_log2FC_curves.csv")
)

p_overall_log2FC <- ggplot(
  overall_log2FC,
  aes(
    x = bin_num,
    y = mean_log2FC,
    color = gene_label,
    fill = gene_label
  )
) +
  geom_ribbon(
    aes(ymin = lower_CI, ymax = upper_CI),
    color = NA,
    alpha = 0.18
  ) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    "nonsig" = "#7852A9",
    "sig" = "#FCAE1E"
  )) +
  scale_fill_manual(values = c(
    "nonsig" = "#7852A9",
    "sig" = "#FCAE1E"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c("-3 Kb", "TSS", "+3 Kb"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = NULL,
    y = expression(log[2]*"((loss + p) / (neutral + p))"),
    color = NULL,
    fill = NULL
  ) +
  ggthemes::geom_rangeframe(color = "black") +
  ggthemes::theme_tufte(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(output_dir, "overall_stabilized_log2FC_sig_vs_nonsig.pdf"),
  p_overall_log2FC,
  width = 6,
  height = 4.5
)

ggsave(
  file.path(output_dir, "overall_stabilized_log2FC_sig_vs_nonsig.png"),
  p_overall_log2FC,
  width = 6,
  height = 4.5,
  dpi = 300
)

# Overall loss and neutral profiles shown in the requested faceted TSS plot.
overall_profile <- df %>%
  group_by(gene_label, bin_num, condition) %>%
  summarise(
    median_value = median(mean, na.rm = TRUE),
    mean_value = mean(mean, na.rm = TRUE),
    n_arms = sum(!is.na(mean)),
    variance_across_arms = var(mean, na.rm = TRUE),
    se_across_arms = sd(mean, na.rm = TRUE) / sqrt(n_arms),
    lower_CI = mean_value - qt(0.975, df = pmax(n_arms - 1, 1)) * se_across_arms,
    upper_CI = mean_value + qt(0.975, df = pmax(n_arms - 1, 1)) * se_across_arms,
    .groups = "drop"
  )

write_csv(
  overall_profile,
  file.path(output_dir, "overall_loss_neutral_profiles.csv")
)

p_overall_profile <- ggplot(
  overall_profile,
  aes(
    x = bin_num,
    y = median_value,
    color = condition,
    fill = condition
  )
) +
  geom_ribbon(
    aes(
      ymin = median_value - se_across_arms,
      ymax = median_value + se_across_arms
    ),
    color = NA,
    alpha = 0.20
  ) +
  geom_smooth(
    method = "loess",
    span = 0.05,
    se = FALSE,
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    "loss" = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    "loss" = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c("-3 Kb", "TSS", "+3 Kb"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = NULL,
    y = "ATAC-seq signal",
    color = NULL,
    fill = NULL
  ) +
  ggthemes::geom_rangeframe(color = "black") +
  ggthemes::theme_tufte(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  facet_wrap(
    ~ gene_label,
    labeller = as_labeller(c(nonsig = "nonsig", sig = "sig"))
  )

ggsave(
  file.path(output_dir, "overall_TSS_profiles_sig_nonsig.pdf"),
  p_overall_profile,
  width = 9,
  height = 4.5
)

ggsave(
  file.path(output_dir, "overall_TSS_profiles_sig_nonsig.png"),
  p_overall_profile,
  width = 9,
  height = 4.5,
  dpi = 300
)

# 2. Difference in stabilized log2FC:
#    log2FC[sig] - log2FC[nonsig]
interaction_df <- difference_df %>%
  select(chr, bin_num, gene_label, difference, variance_difference) %>%
  pivot_wider(
    names_from = gene_label,
    values_from = c(difference, variance_difference)
  ) %>%
  mutate(
    sig_log2FC = difference_sig,
    nonsig_log2FC = difference_nonsig,
    DID = difference_sig - difference_nonsig,
    log2FC_contrast = DID,
    variance_DID = variance_difference_sig + variance_difference_nonsig,
    se_DID = sqrt(variance_DID),
    z_value = DID / se_DID,
    p_pointwise = 2 * pnorm(-abs(z_value)),
    lower_CI = DID - 1.96 * se_DID,
    upper_CI = DID + 1.96 * se_DID
  ) %>%
  group_by(chr) %>%
  mutate(FDR_within_chr = p.adjust(p_pointwise, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_all_bins = p.adjust(p_pointwise, method = "BH"))

write_csv(
  interaction_df,
  file.path(output_dir, "difference_in_differences_by_bin.csv")
)

# Overall sig-minus-nonsig log2FC curve across chromosome arms.
# Negative values mean a stronger reduction in significant paralogs.
# Positive values mean a stronger increase (or weaker reduction) in sig paralogs.
overall_interaction <- interaction_df %>%
  group_by(bin_num) %>%
  summarise(
    mean_DID = mean(DID, na.rm = TRUE),
    median_DID = median(DID, na.rm = TRUE),
    n_arms = sum(!is.na(DID)),
    variance_across_arms = var(DID, na.rm = TRUE),
    se_mean_DID = sqrt(variance_across_arms / n_arms),
    lower_CI = mean_DID - qt(0.975, df = pmax(n_arms - 1, 1)) * se_mean_DID,
    upper_CI = mean_DID + qt(0.975, df = pmax(n_arms - 1, 1)) * se_mean_DID,
    t_value = mean_DID / se_mean_DID,
    p_pointwise = 2 * pt(-abs(t_value), df = pmax(n_arms - 1, 1)),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_pointwise, method = "BH"))

write_csv(
  overall_interaction,
  file.path(output_dir, "overall_sig_minus_nonsig_log2FC.csv")
)

# Global functional sign-flip test. All bins stay together during permutation,
# so correlation between neighboring bins is preserved. The global statistic is
# the mean squared value of the overall curve.
functional_signflip_test <- function(input_df, value_column, n_perm = 10000,
                                     seed = 20260916) {
  value_column <- rlang::ensym(value_column)

  curve_matrix <- input_df %>%
    select(chr, bin_num, value = !!value_column) %>%
    pivot_wider(names_from = bin_num, values_from = value) %>%
    arrange(chr)

  arm_names <- curve_matrix$chr
  mat <- as.matrix(select(curve_matrix, -chr))
  observed_curve <- colMeans(mat, na.rm = TRUE)
  observed_statistic <- mean(observed_curve^2, na.rm = TRUE)

  set.seed(seed)
  permuted_statistics <- replicate(n_perm, {
    signs <- sample(c(-1, 1), nrow(mat), replace = TRUE)
    permuted_curve <- colMeans(mat * signs, na.rm = TRUE)
    mean(permuted_curve^2, na.rm = TRUE)
  })

  tibble(
    n_arms = length(arm_names),
    n_permutations = n_perm,
    observed_global_statistic = observed_statistic,
    global_p = (sum(permuted_statistics >= observed_statistic) + 1) /
      (n_perm + 1)
  )
}

# Global tests for the two loss-neutral curves and their interaction.
global_nonsig <- difference_df %>%
  filter(gene_label == "nonsig") %>%
  functional_signflip_test(difference) %>%
  mutate(comparison = "stabilized log2FC loss versus neutral: nonsig")

global_sig <- difference_df %>%
  filter(gene_label == "sig") %>%
  functional_signflip_test(difference) %>%
  mutate(comparison = "stabilized log2FC loss versus neutral: sig")

global_interaction <- interaction_df %>%
  functional_signflip_test(DID) %>%
  mutate(comparison = "stabilized log2FC[sig] - stabilized log2FC[nonsig]")

overall_global_tests <- bind_rows(
  global_nonsig,
  global_sig,
  global_interaction
) %>%
  select(comparison, everything()) %>%
  mutate(FDR = p.adjust(global_p, method = "BH"))

write_csv(
  overall_global_tests,
  file.path(output_dir, "overall_global_functional_tests.csv")
)

overall_effect_size <- overall_interaction %>%
  arrange(bin_num) %>%
  summarise(
    overall_mean_DID = mean(mean_DID, na.rm = TRUE),
    mean_absolute_DID = mean(abs(mean_DID), na.rm = TRUE),
    RMS_DID = sqrt(mean(mean_DID^2, na.rm = TRUE)),
    maximum_absolute_DID = max(abs(mean_DID), na.rm = TRUE),
    bin_of_maximum = bin_num[which.max(abs(mean_DID))],
    mean_variance_across_arms = mean(variance_across_arms, na.rm = TRUE),
    signed_AUC = sum(
      (head(mean_DID, -1) + tail(mean_DID, -1)) / 2 * diff(bin_num),
      na.rm = TRUE
    ),
    absolute_AUC = sum(
      (head(abs(mean_DID), -1) + tail(abs(mean_DID), -1)) / 2 * diff(bin_num),
      na.rm = TRUE
    )
  )

write_csv(
  overall_effect_size,
  file.path(output_dir, "overall_interaction_effect_size.csv")
)

p_overall_did <- ggplot(
  overall_interaction,
  aes(x = bin_num, y = mean_DID)
) +
  geom_ribbon(
    aes(ymin = lower_CI, ymax = upper_CI),
    fill = "#7852A9",
    alpha = 0.20
  ) +
  geom_line(color = "#7852A9", linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = tss_pos, linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c("-3 Kb", "TSS", "+3 Kb"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = NULL,
    y = "Sig - nonsig stabilized log2FC"
  ) +
  ggthemes::geom_rangeframe(color = "black") +
  ggthemes::theme_tufte(base_size = 12)

ggsave(
  file.path(output_dir, "overall_sig_minus_nonsig_log2FC.pdf"),
  p_overall_did,
  width = 6,
  height = 4.5
)

ggsave(
  file.path(output_dir, "overall_sig_minus_nonsig_log2FC.png"),
  p_overall_did,
  width = 6,
  height = 4.5,
  dpi = 300
)

# 3. Global GAM interaction test for each chromosome.
# Null: sig and nonsig have the same stabilized log2FC curve shape.
# Full: sig has an additional smooth deviation from the nonsig curve.
run_gam_test <- function(dat) {
  model_df <- dat %>%
    mutate(
      gene_label = factor(gene_label, levels = c("nonsig", "sig")),
      is_sig = as.numeric(gene_label == "sig"),
      weight = 1 / pmax(variance_difference, 1e-6)
    ) %>%
    filter(
      is.finite(difference),
      is.finite(bin_num),
      is.finite(weight)
    )

  fit_null <- gam(
    difference ~ gene_label + s(bin_num, k = 20, bs = "cs"),
    data = model_df,
    weights = weight,
    method = "ML"
  )

  fit_full <- gam(
    difference ~ gene_label +
      s(bin_num, k = 20, bs = "cs") +
      s(bin_num, by = is_sig, k = 20, bs = "cs"),
    data = model_df,
    weights = weight,
    method = "ML"
  )

  model_comparison <- anova(fit_null, fit_full, test = "Chisq")
  p_column <- grep("^Pr", names(model_comparison), value = TRUE)[1]
  global_p <- if (length(p_column) == 1) model_comparison[[p_column]][2] else NA_real_
  full_summary <- summary(fit_full)

  tibble(
    global_curve_p = global_p,
    AIC_null = AIC(fit_null),
    AIC_full = AIC(fit_full),
    deviance_explained = full_summary$dev.expl,
    residual_variance = full_summary$scale,
    n_bins_used = nrow(model_df)
  )
}

gam_results <- difference_df %>%
  group_by(chr) %>%
  group_modify(~ run_gam_test(.x)) %>%
  ungroup() %>%
  mutate(global_FDR = p.adjust(global_curve_p, method = "BH"))

# 4. Effect sizes for the sig-minus-nonsig stabilized-log2FC contrast.
# Since bins are equally spaced, AUC is reported in log2FC x bin units.
effect_sizes <- interaction_df %>%
  arrange(chr, bin_num) %>%
  group_by(chr) %>%
  summarise(
    mean_DID = mean(DID, na.rm = TRUE),
    mean_absolute_DID = mean(abs(DID), na.rm = TRUE),
    RMS_DID = sqrt(mean(DID^2, na.rm = TRUE)),
    maximum_absolute_DID = max(abs(DID), na.rm = TRUE),
    bin_of_maximum = bin_num[which.max(abs(DID))],
    signed_AUC = sum(
      (head(DID, -1) + tail(DID, -1)) / 2 * diff(bin_num),
      na.rm = TRUE
    ),
    absolute_AUC = sum(
      (head(abs(DID), -1) + tail(abs(DID), -1)) / 2 * diff(bin_num),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

summary_results <- gam_results %>%
  left_join(effect_sizes, by = "chr") %>%
  arrange(global_FDR)

write_csv(
  summary_results,
  file.path(output_dir, "global_curve_tests_and_effect_sizes.csv")
)

# 5. Plot chromosome-specific sig-minus-nonsig stabilized-log2FC curves.
p_did <- ggplot(
  interaction_df,
  aes(x = bin_num, y = DID)
) +
  geom_ribbon(
    aes(ymin = lower_CI, ymax = upper_CI),
    fill = "#7852A9",
    alpha = 0.20
  ) +
  geom_line(color = "#7852A9", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = tss_pos, linetype = "dashed", color = "black") +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c("-3 Kb", "TSS", "+3 Kb"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  facet_wrap(~ chr, scales = "free_y", ncol = 4) +
  labs(
    x = NULL,
    y = "Sig - nonsig stabilized log2FC"
  ) +
  theme_classic() +
  theme(strip.background = element_blank())

ggsave(
  file.path(output_dir, "chromosome_stabilized_log2FC_contrasts.pdf"),
  p_did,
  width = 12,
  height = 11
)

ggsave(
  file.path(output_dir, "chromosome_stabilized_log2FC_contrasts.png"),
  p_did,
  width = 12,
  height = 11,
  dpi = 300
)

# 6. Summary plot: global effect magnitude and adjusted P value.
p_summary <- summary_results %>%
  mutate(chr = reorder(chr, mean_absolute_DID)) %>%
  ggplot(aes(x = chr, y = mean_absolute_DID)) +
  geom_col(aes(fill = global_FDR < 0.05), width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = c(`TRUE` = "#FCAE1E", `FALSE` = "grey70"),
    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05")
  ) +
  labs(
    x = NULL,
    y = "Mean absolute stabilized-log2FC difference",
    fill = NULL
  ) +
  theme_classic()

ggsave(
  file.path(output_dir, "global_effect_size_summary.pdf"),
  p_summary,
  width = 7,
  height = 6
)

print(summary_results)


####### look at individual chromosome arms 
## compare the sig and non-sig paralog , the signal between (loss and neutral)

interaction_df %>% 
  dplyr::filter(chr %in% c("chr8p", "chr11q", 
                           "chr17","chr5q" )) %>% 
  dplyr::mutate(chr = factor(chr,
                             levels = c("chr8p", "chr5q" ,
                                        "chr11q", "chr17"))) %>% 
  ggplot(
 
  aes(x = bin_num, y = DID)
  ) +
  geom_ribbon(
    aes(ymin = lower_CI, ymax = upper_CI),
    fill = "#7852A9",
    alpha = 0.20
  ) +
  geom_line(color = "#7852A9", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = tss_pos, linetype = "dashed", color = "black") +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c("-3 Kb", "TSS", "+3 Kb"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  facet_wrap(~ chr, scales = "free_y", ncol = 4) +
  labs(
    x = NULL,
    y = "Sig - nonsig stabilized log2FC"
  ) +
  geom_rangeframe() +
  theme_tufte() +
  theme(strip.background = element_blank())
# Statistical note:
# Pointwise variances assume the four means are independent. The GAM P values
# treat summarized bins as observations and therefore do not fully model the
# correlation between neighboring bins. For definitive inference, repeat this
# analysis using gene-level or biological-sample-level curves and permute labels
# at that level.
