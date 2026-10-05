# Project: H3K27ac HiChIP Peaks Across Recurrent Chromosome-Loss Events in Cancer
# Date: Sep.10.2026
# Author: Yi Bei

# load the library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# laod the data from the HiChIP
# The original paper is from "Three-dimensional genome landscape of primary human cancers, 2025"

# Load the data
# Metasheet data for the HiChIP
# get all the sample info 

sample_meta_HiChIP = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/HiChIP dataset/TCGA_HiChIP_metadata.txt")

# TCGA sample aneuploidy
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
sample_folder = list.files("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/HiChIP dataset/H3K27ac_all_chromosome",full.names = T)



for (folder in sample_folder) {
  
  ### get the files
  files <- list.files(
    path = folder,
    pattern = "\\.tab$",
    full.names = TRUE)
  
  # get the chromosome location
  chr <- str_extract(folder, "(?:[1-9]|1[0-9]|2[0-2]|X|Y)[pq]")
  
  if(is.na(chr)){
    chr = "17p"
  }
  
  
  print(chr)
  # get the loss sample
  
  # get loss chromosome 
  chr_loss_df = TACG_df[TACG_df[,chr] == -1 &TACG_df$Genome_doublings == 0 ,]
  
  chr_loss_sample = chr_loss_df$Sample
  chr_loss_sample = chr_loss_sample[!is.na(chr_loss_sample)]
  
  # to get the ATACseq_sample
  sample_meta_HiChIP %>% 
    dplyr::mutate(sample_id = str_sub(aliquot_submitter_id, 1,15)) %>% 
    dplyr::filter(sample_id %in% chr_loss_sample) %>% 
    dplyr::mutate(bam_prefix = sub(
      "(_H3K27ac).*$",
      "\\1",
      basename(file_name)
    )) %>% 
    dplyr::select(bam_prefix) %>% 
    distinct(bam_prefix) %>% 
    unlist(use.names = F) -> chr_loss_hiChip_sample
  
  if(length(chr_loss_hiChip_sample) >= 5){
    # get the loss files
    files_chr_loss_selected <- files[basename(files) %in% paste(chr_loss_hiChip_sample, "_1D-signal.tab", sep = "")]
    
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
    sample_meta_HiChIP %>% 
      dplyr::mutate(sample_id = str_sub(aliquot_submitter_id, 1,15)) %>% 
      dplyr::filter(sample_id %in% chr_neutral_sample) %>% 
      dplyr::mutate(bam_prefix = sub(
        "(_H3K27ac).*$",
        "\\1",
        basename(file_name)
      )) %>% 
      dplyr::select(bam_prefix) %>% 
      distinct(bam_prefix) %>% 
      unlist(use.names = F) -> chr_neutral_hiChip_sample
    
    # get the loss files
    files_chr_neutral_selected <- files[basename(files) %in% paste(chr_neutral_hiChip_sample, "_1D-signal.tab", sep = "")]
    
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
    
    
    
    # add the gene label
    df_cmb$gene_label = sub("^.*_(.*?)\\.bed$", "\\1", df_cmb$group)
    
    # add chromosome information
    df_cmb$chr = sub("_.*$", "", df_cmb$group)
    # check the tbl
    head(df_cmb)
    
    # save the file into big formation
    if(!exists("comb_df")){
      comb_df = df_cmb
    }else{
      comb_df = rbind(comb_df, df_cmb)
    }
    
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
# dplyr::filter(gene_label == "chr17")


#The missing ribbon occurs because sd() returns NA when a bin has only one observation.
#Also, completely absent bins do not appear in plot_df. 
#Complete all bins first, then interpolate missing median and SE values within each gene_label × condition group.
interpolate_na <- function(x, y) {
  valid <- is.finite(y)
  
  if (sum(valid) >= 2) {
    approx(
      x = x[valid],
      y = y[valid],
      xout = x,
      rule = 2
    )$y
  } else if (sum(valid) == 1) {
    rep(y[valid][1], length(y))
  } else {
    rep(NA_real_, length(y))
  }
}

plot_df <- comb_df %>%
  dplyr::group_by(gene_label, bin, condition) %>%
  dplyr::summarise(
    median_value = median(mean, na.rm = TRUE),
    n_files = sum(!is.na(mean)),
    
    # SE cannot be calculated from only one sample
    se_value = dplyr::if_else(
      n_files >= 2,
      sd(mean, na.rm = TRUE) / sqrt(n_files),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    median_value = dplyr::na_if(median_value, NaN)
  ) %>%
  
  # Add bins that are completely absent from the original data
  dplyr::group_by(gene_label, condition) %>%
  tidyr::complete(bin = seq_len(total_bins)) %>%
  dplyr::arrange(bin, .by_group = TRUE) %>%
  
  # Interpolate missing values within each curve
  dplyr::mutate(
    median_value = interpolate_na(bin, median_value),
    se_value = interpolate_na(bin, se_value),
    ribbon_lower = pmax(0, median_value - se_value),
    ribbon_upper = median_value + se_value
  ) %>%
  dplyr::ungroup()



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
      ymin = ribbon_lower,
      ymax = ribbon_upper,
      group = condition
    ),
    color = NA,
    alpha = 0.20
  ) +
  geom_smooth(
    aes(group = condition),
    method = "loess",
    span = 0.05,
    se = FALSE,
    linewidth = 1.1
  ) +
  geom_vline(
    xintercept = tss_pos,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "loss" = "#FCAE1E",
      "neutral" = "#7852A9"
    )
  ) +
  scale_fill_manual(
    values = c(
      "loss" = "#FCAE1E",
      "neutral" = "#7852A9"
    )
  ) +
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
  geom_rangeframe(color = "black") +
  theme_tufte(base_size = 12) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  ) +
  facet_wrap(~gene_label)

# 800, 300

head(comb_df)
#writexl::write_xlsx(comb_df, "RChIP_H3K27ac_TSS.xlsx")
