# Project: Plot alll the other chromosome loss events ATACseq
# Date: Sep_09_2026
# Author: Yi

# library
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

##data
sample_atac = read.delim("../TCGA_identifier_mapping.txt")
TACG_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Aneuploidy patient datasets/mmc2.xlsx", skip = 1)


# function 
# ----------------------------
# function to read one matrix
# ----------------------------
read_profile_tss <- function(file) {
  
  x <- readLines(file)
  
  # line 1: group names and group sizes
  h1 <- sub("^#", "", x[1])
  ginfo <- strsplit(h1, "\t")[[1]]
  
  gnames <- sub(":(.*)$", "", ginfo)
  gsizes <- as.numeric(sub("^.*:", "", ginfo))
  
  # line 2: parameters
  h2 <- sub("^#", "", x[2])
  
  upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
  downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
  bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))
  
  # numeric matrix starts after first 3 lines
  num_lines <- x[-c(1, 2, 3)]
  
  mat <- do.call(
    rbind,
    lapply(num_lines, function(z) as.numeric(strsplit(z, "\t")[[1]]))
  )
  
  mat <- as.matrix(mat)
  
  group <- rep(gnames, gsizes)
  
  out <- do.call(
    rbind,
    lapply(gnames, function(g) {
      idx <- which(group == g)
      
      data.frame(
        group = g,
        bin = seq_len(ncol(mat)),
        signal = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
      )
    })
  )
  
  out$sample <- basename(file)
  out$upstream <- upstream
  out$downstream <- downstream
  out$bin_size <- bin_size
  
  out
}
### get the folder sample ###
sample_folder = list.files("../Outputs_data",full.names = T)

for (folder in sample_folder) {
  
  ### get the files
  files <- list.files(
    path = folder,
    pattern = "\\.tab$",
    full.names = TRUE)
  
  # get the chromosome location
  chr = str_sub(folder, start = 27, end = nchar(folder))
  
  # get the loss sample
  
  # get loss chromosome 
  chr_loss_df = TACG_df[TACG_df[,chr] == -1 &TACG_df$Genome_doublings == 0 ,]
  
  chr_loss_sample = chr_loss_df$Sample
  chr_loss_sample = chr_loss_sample[!is.na(chr_loss_sample)]
  
  # to get the ATACseq_sample
  sample_atac %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_loss_sample) %>% 
    dplyr::mutate(bam_prefix = gsub("-", "_", bam_prefix)) %>% 
    dplyr::select(bam_prefix) %>% 
    unlist(use.names = F) -> chr_loss_atac_sample
  
  # get the loss files
  files_chr_loss_selected <- files[basename(files) %in% paste(chr_loss_atac_sample, ".tab", sep = "")]
  
  # read the loss files 
  chr_loss_profiles <- bind_rows(lapply(files_chr_loss_selected, read_profile_tss))
  
  # ----------------------------
  # average across loss files
  # ----------------------------
  summary_loss_df <- chr_loss_profiles %>%
    group_by(group, bin) %>%
    summarise(
      mean = mean(signal, na.rm = TRUE),
      sd   = sd(signal, na.rm = TRUE),
      n    = n(),
      se   = sd / sqrt(n),
      .groups = "drop"
    )
  
  
  # get the neutral files 
  # get loss chromosome 
  chr_neutral_df = TACG_df[TACG_df[,chr] == 0 &TACG_df$Genome_doublings == 0 ,]
  
  chr_neutral_sample = chr_neutral_df$Sample
  chr_neutral_sample = chr_neutral_sample[!is.na(chr_neutral_sample)]
  
  # to get the ATACseq_sample
  sample_atac %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_neutral_sample) %>% 
    dplyr::mutate(bam_prefix = gsub("-", "_", bam_prefix)) %>% 
    dplyr::select(bam_prefix) %>% 
    unlist(use.names = F) -> chr_neutral_atac_sample
  
  # get the loss files
  files_chr_neutral_selected <- files[basename(files) %in% paste(chr_neutral_atac_sample, ".tab", sep = "")]
  
  # read the loss files 
  chr_neutral_profiles <- bind_rows(lapply(files_chr_neutral_selected, read_profile_tss))
  
  
  # get all the neutral files 
  summary_neutral_df <- chr_neutral_profiles %>%
    group_by(group, bin) %>%
    summarise(
      mean = mean(signal, na.rm = TRUE),
      sd   = sd(signal, na.rm = TRUE),
      n    = n(),
      se   = sd / sqrt(n),
      .groups = "drop"
    )
  
  
  # normalise the copy number
  summary_neutral_df %>% 
    dplyr::mutate(mean = mean/ 2) -> summary_neutral_df_mdf
  
  summary_neutral_df_mdf$condition = "neutral"
  summary_loss_df$condition = "loss"
  
  # combine the table 
  df_cmb = rbind(summary_neutral_df_mdf,summary_loss_df)
  
  # 8p file is different 
  if(chr == "8p"){
    df_cmb$chr = "chr8p"
    df_cmb %>% 
      mutate(gene_label = ifelse(group == "PPP2CB.bed", "sig", "nonsig")) ->df_cmb
  }else{
    # add the gene label
    df_cmb$gene_label = sub("^.*_(.*?)\\.bed$", "\\1", df_cmb$group)
    
    # add chromosome information
    df_cmb$chr = sub("_.*$", "", df_cmb$group)
    # check the tbl
    head(df_cmb)
  }
 
  
  # save the file into big formation
  if(!exists("comb_df")){
    comb_df = df_cmb
  }else{
    comb_df = rbind(comb_df, df_cmb)
  }
  
}


h2 <- sub("^#", "", readLines(files_chr_loss_selected[1], n = 2)[2])

upstream   <- as.numeric(sub(".*upstream:(\\d+).*", "\\1", h2))
downstream <- as.numeric(sub(".*downstream:(\\d+).*", "\\1", h2))
bin_size   <- as.numeric(sub(".*bin size:(\\d+).*", "\\1", h2))

up_bins   <- upstream / bin_size
down_bins <- downstream / bin_size
total_bins <- up_bins + down_bins

tss_pos <- up_bins + 1

#comb_df_sig <- comb_df %>%
 # dplyr::filter(gene_label == "sig")

plot_df <- comb_df %>%
  #dplyr::filter(chr == "chr8p")
  dplyr::group_by(gene_label, bin, condition) %>%
  dplyr::summarise(
    median_value = median(mean, na.rm = TRUE),
    n_files      = sum(!is.na(mean)),
    se_value     = sd(mean, na.rm = TRUE) / sqrt(n_files),
    .groups = "drop"
  )


# plot the figure
library(ggthemes)
ggplot(
  plot_df,
  aes(
    x = bin,
    y = median_value,
    color = condition,
    fill = condition
  )
) +
  geom_ribbon(
    aes(
      ymin = median_value - se_value,
      ymax = median_value + se_value
    ),
    color = NA,
    alpha = 0.20
  ) +
  geom_smooth(
    method = "loess",
    span = 0.05,
    se = FALSE,
    linewidth = 1.1
  )+
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    "loss"    = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    "loss"    = "#FCAE1E",
    "neutral" = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  #coord_cartesian(ylim = c(0, 25)) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  geom_rangeframe(color = "black") +
  theme_tufte(base_size = 12)+
  #theme_classic(base_size = 18) +
  theme(
  #  legend.position = c(0.5, 0.85),
     #axis.line = element_line(color =  "black"),
  #  axis.ticks = element_line(linewidth = 0.8),
      plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  facet_wrap(~gene_label)


# 800, 300

####### look at the inidvidual chromosome arm loss ##########
# Chr8p loss#
plot_df_8p <- comb_df %>%
  filter(chr == "chr8p", is.finite(bin), is.finite(mean)) %>%
  group_by(gene_label, condition) %>%
  group_modify(~ {
    df <- arrange(.x, bin)
    
    fit <- loess(
      mean ~ bin,
      data = df,
      span = 0.05,
      na.action = na.exclude
    )
    
    df %>%
      mutate(
        mean_smooth = as.numeric(
          predict(fit, newdata = data.frame(bin = bin))
        ),
        lower = mean_smooth - 1.5 * se,
        upper = mean_smooth + 1.5 * se
      )
  }) %>%
  ungroup()

ggplot(
  plot_df_8p,
  aes(
    x = bin,
    y = mean_smooth,
    color = condition,
    fill = condition,
    group = condition
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    color = NA,
    alpha = 0.20
  ) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  ggthemes::geom_rangeframe(color = "black") +
  ggthemes::theme_tufte(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  facet_wrap(~gene_label)


########### 5q ###########
plot_df_5q <- comb_df %>%
  filter(chr == "chr5q", is.finite(bin), is.finite(mean)) %>%
  group_by(gene_label, condition) %>%
  group_modify(~ {
    df <- arrange(.x, bin)
    
    fit <- loess(
      mean ~ bin,
      data = df,
      span = 0.05,
      na.action = na.exclude
    )
    
    df %>%
      mutate(
        mean_smooth = as.numeric(
          predict(fit, newdata = data.frame(bin = bin))
        ),
        lower = mean_smooth - 1.5 * se,
        upper = mean_smooth + 1.5 * se
      )
  }) %>%
  ungroup()

ggplot(
  plot_df_5q,
  aes(
    x = bin,
    y = mean_smooth,
    color = condition,
    fill = condition,
    group = condition
  )
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    color = NA,
    alpha = 0.20
  ) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_fill_manual(values = c(
    loss = "#FCAE1E",
    neutral = "#7852A9"
  )) +
  scale_x_continuous(
    breaks = c(1, tss_pos, total_bins),
    labels = c(
      paste0("-", upstream / 1000, " Kb"),
      "TSS",
      paste0("+", downstream / 1000, " Kb")
    )
  ) +
  labs(
    x = NULL,
    y = "Signal",
    color = NULL,
    fill = NULL
  ) +
  ggthemes::geom_rangeframe(color = "black") +
  ggthemes::theme_tufte(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  facet_wrap(~gene_label)
