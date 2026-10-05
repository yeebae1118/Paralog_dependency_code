# Project: paralog pair protein location frequency 
# Author: Yi Bei
# Date: May.27.2026

# Library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# dataset
localisation_protein_df = read.delim("protein_localisation_dataset/subcellular_location.tsv")

# significant paralog
sig_paralog = readxl::read_xlsx("protein_localisation_dataset/sig_37_paralog.xlsx")

# non-significant paralog 
non_sig_paralog = readxl::read_xlsx("protein_localisation_dataset/non_sig_paralog.xlsx")

# first we have to know how many pair of can not find in the protein dataset in the protein atlas 
# Set rule: 
# either of paralog or neither of paralog gene in the dataset , we will exclude the analysis afterwards

# function 
# index here as input data
OVERLAPPED_FREQ = function(data, index_length){
  
  # generate vector
  frequency =  c()
  # get the input data
  for (i in c(1:index_length)) {
    # get the data by row
    df_tpm = data[i,]
    
    # get the paralog_1
    para_1 = unlist(df_tpm[,1],use.names = F)
    
    # get the paralog_2
    para_2 = unlist(df_tpm[,2],use.names = F)
    
    # look in to the location data sheet and get the location of protein (approved or supported)
    # both conditions are included into our analysis 
    
    # location of paralog1
    localisation_protein_df %>% 
      dplyr::filter(Gene.name == para_1) -> loc_flt_para_1
    
    # get the main positionand additional position 
    position_para_1 = c(loc_flt_para_1$Main.location,loc_flt_para_1$Additional.location)
    position_para_1_splt = unlist(strsplit(position_para_1, ";"))
    
    # location of paralog2
    localisation_protein_df %>% 
      dplyr::filter(Gene.name == para_2) -> loc_flt_para_2
    
    # get the main positionand additional position 
    position_para_2 = c(loc_flt_para_2$Main.location,loc_flt_para_2$Additional.location)
    position_para_2_splt = unlist(strsplit(position_para_2, ";"))
    
    # check for the intersection 
    overlapped_region = intersect(position_para_1_splt, position_para_2_splt)
    
    # calculate the mean frequency
    mean_freq = (length(overlapped_region)/length(position_para_1_splt) + 
                   length(overlapped_region)/length(position_para_2_splt)) / 2
    
    frequency = c(frequency, mean_freq)  
  }
  
  return(frequency)
}

# first check for the significant paralog pairs
sig_count = 0
sig_position = c()
for (i in c(1: dim(sig_paralog)[1])) {
  
  # get the df
  df_tmp = sig_paralog[i,]
  # get the paralog gene
  para_1 = unlist(df_tmp[,1],use.names = F)
  
  # get the paralog partner 
  para_2 = unlist(df_tmp[,7],use.names = F)
  
  # check if both of gene in the database of protein localation dataset
  # check for the parlaog 1
  sum_1 = sum(para_1 %in% localisation_protein_df$Gene.name)
  
  # check for the paralog 2
  sum_2 = sum(para_2 %in% localisation_protein_df$Gene.name)
  
  # if sum1 + sum2 = 2 means both sample are presenting in the dataset 
  if(sum_1 + sum_2 == 2){
    sig_count = sig_count + 1
    
    # add into the position list 
    sig_position = c(sig_position, i)
  }

   
}

# check for the non-significant 
non_sig_count = 0
non_sig_position = c()
for (i in c(1: dim(non_sig_paralog)[1])) {
  
  # get the df
  df_tmp = non_sig_paralog[i,]
  # get the paralog gene
  para_1 = unlist(df_tmp[,3],use.names = F)
  
  # get the paralog partner 
  para_2 = unlist(df_tmp[,4],use.names = F)
  
  # check if both of gene in the database of protein localation dataset
  # check for the parlaog 1
  sum_1 = sum(para_1 %in% localisation_protein_df$Gene.name)
  
  # check for the paralog 2
  sum_2 = sum(para_2 %in% localisation_protein_df$Gene.name)
  
  # if sum1 + sum2 = 2 means both sample are presenting in the dataset 
  if(sum_1 + sum_2 == 2){
    non_sig_count = non_sig_count + 1
    
    # add into the position list 
    non_sig_position = c(non_sig_position, i)
  }
  
  
}

# plot the presenting of sig and non-sig paralog pair presenting in the dataset 
# 25/37
# 13062/28103
summary_df = tibble::tibble(freq = c(sig_count/dim(sig_paralog)[1], non_sig_count/dim(non_sig_paralog)[1]),
                            annotation = c("significant paralog", "non-significant paralog"))

# plot the statistic 
summary_df %>% 
  ggplot(aes(x = annotation, y = freq)) +
  geom_col(aes(color = annotation, 
               fill = annotation), width = 0.5, alpha = 0.7) +
  theme_tufte() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# check for the position that is significant 
sig_flt_df = sig_paralog[sig_position,c(1,7)]


# to check for the overlapping frequency in the significant pairs
freq_sig = OVERLAPPED_FREQ(sig_flt_df, dim(sig_flt_df)[1]) 
freq_df_sig = data.frame(freq_sig)
freq_df_sig$annotation = "sig"
colnames(freq_df_sig)[1] = "freq"
# to check for the non-sig freq
nonsig_flt_df = non_sig_paralog[non_sig_position,c(3,4)]

# to check for the overlapping frequency in the non_significant pairs
freq_non_sig = OVERLAPPED_FREQ(nonsig_flt_df, dim(nonsig_flt_df)[1]) 

freq_df_nonsig = data.frame(freq_non_sig)
freq_df_nonsig$annotation = "non-sig"
colnames(freq_df_nonsig)[1] = "freq"

# combine the data
freq_df_cmb = rbind(freq_df_sig, freq_df_nonsig)

colors = c("#EF8536","#3A76AF")
# plot the figure
freq_df_cmb %>% 
  mutate(annotation = factor(annotation, levels = c("sig", "non-sig"))) %>% 
  ggplot(aes(x = annotation, y = freq)) +
  geom_boxplot(width = 0.25, 
               aes(color = annotation, fill = annotation),
               alpha = 0.5,
               outliers = F) +
  #geom_jitter( aes(color = annotation),position=position_jitter(0.12),size = 1,alpha = 0.15) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[1], "#7852A9"))+
  scale_fill_manual(values = c(colors[1], "#7852A9"))
  
####### t.test ##########
# p-value = 0.01451
t.test(freq_df_cmb[freq_df_cmb$annotation == "sig",]$freq, 
       freq_df_cmb[freq_df_cmb$annotation == "non-sig",]$freq)
