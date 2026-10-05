# Project: CYCLOPS are singleton 
# Date: Dec.8th
# Author: Yi

# library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# colors
colors = c("#EF8536","#3A76AF")
# dataset
# load hg38 genomes
library(biomaRt)

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

mart <- useMart(biomart="ensembl",
                dataset="hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position",
                "end_position", "strand", "gene_biotype", "description","band")

# Fetch the gene information
gene_info <- getBM(attributes = attributes, mart = mart,useCache = TRUE)

library(MetBrewer)
# check the gene annotation and distribution
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22,"Y","X")) %>% 
  dplyr::filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% 
  group_by(gene_biotype) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  dplyr::mutate(gene_biotype = factor(gene_biotype, levels = gene_biotype))->p
p %>% 
  ggplot(aes(x =gene_biotype, y = count)) +
  geom_histogram(stat = "identity", width = 0.7,aes(color = gene_biotype, fill = gene_biotype), alpha = 0.8) +
  geom_rangeframe() +
  theme_tufte()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "none") +
  scale_fill_manual(values  = met.brewer("Hokusai3", n = 33 ,direction = 1)) +
  scale_color_manual(values = met.brewer("Hokusai3" ,n = 33,direction = 1))

# include all the genes from the 
# check how many genes with paralogs
# focus on the with protein-coding
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22,"Y","X")) %>% 
  dplyr::filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% 
  dplyr::filter(gene_biotype == "protein_coding")-> all_chr_gene


# gene with parlogus 
paralogs <- getBM(attributes = c("hsapiens_paralog_associated_gene_name",
                                 "ensembl_gene_id",
                                 "hsapiens_paralog_ensembl_gene",
                                 "hsapiens_paralog_perc_id"),
                  filters = "ensembl_gene_id", 
                  values =all_chr_gene$ensembl_gene_id, 
                  mart = mart)


# integrate the genes_names with paralogs data
# filter out the no identical score pair
# get the gene with pairs
all_chr_gene %>% 
  left_join(paralogs, by = "ensembl_gene_id") %>% 
  dplyr::filter(!is.na(hsapiens_paralog_perc_id)) %>% 
  dplyr::filter(hsapiens_paralog_associated_gene_name != "") %>% 
  dplyr::mutate(para_pair_for = paste(hgnc_symbol,hsapiens_paralog_associated_gene_name,sep = "_" )) %>% 
  dplyr::mutate(para_pair_rev = paste(hsapiens_paralog_associated_gene_name,hgnc_symbol,sep = "_" ))-> paralog_gene_pair

# check the head of paralogs pair
View(head(paralog_gene_pair))

# there is an issue that the paralogs have the forward and reverse pair
# this goona cause the duplicate he count of genes
# we need to avoid the duplicate count of paralogs 
# identical score is different calculate in the same pair of paralogs 
# we goona gonna calculate the mean of identical score

for_para = paralog_gene_pair$para_pair_for
rev_para = paralog_gene_pair$para_pair_rev

# read the paralog pair both side
reformed_para_cmb = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/paralog_pair_both_side.csv")
reformed_para_cmb = reformed_para_cmb[,-1]

## To get the unique paralog pairs
# remove the duplicates
# index_rev corresponding on paralog_gene_pair
reformed_para_cmb %>% 
  dplyr::filter(index_for < index_rev | is.na(index_rev) ) %>% 
  distinct(para_name, .keep_all = T)-> reformed_para_cmb_flt # distinct the duplication # no duplicates paralog pair

# define the paralgs by identical score
# calcluate the famaily size
reformed_para_cmb %>% # paralog gene family
  distinct(para_name, .keep_all = T) %>% 
  dplyr::filter(mean_identical_score >=20) %>% 
  group_by(hgnc_symbol) %>% 
  summarise(family_size = n()) -> family_size_sum

# filter the paralogs
reformed_para_cmb_flt %>% 
  dplyr::filter(mean_identical_score >= 20) %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  left_join(family_size_sum, by = "hgnc_symbol") %>% 
  dplyr::filter(family_size < 20) -> paralog_pair_flt


# get the CYClOPS gene
# read the significant target 
sig_target = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")

sig_target %>% 
  dplyr::filter(status == "CYCLOPS") -> df_CYCLOPS

#
paralog_num = c()
for (gene in df_CYCLOPS$dep_gene) {
  if(!(gene %in% family_size_sum$hgnc_symbol)){
    family_size = 0
    # add the family size
    paralog_num = c(paralog_num, family_size)
  }else{
    family_size_sum %>% 
      dplyr::filter(hgnc_symbol == gene) %>% 
      distinct(hgnc_symbol, .keep_all = T)-> tpm_df
    family_size = tpm_df$family_size
    # add the family size
    paralog_num = c(paralog_num, family_size)
    
    
  }
}

CYCLOPS_family_size = data.frame(df_CYCLOPS$dep_gene, paralog_num)

CYCLOPS_family_size %>% 
  dplyr::mutate(paralog_num = paralog_num +1) -> CYCLOPS_family_size_mdf
CYCLOPS_family_size_mdf %>% 
  ggplot(aes(x = paralog_num)) +
  geom_bar(color = "#5271A8", fill ="#5271A8", alpha = 0.6,width = 0.75) +
  geom_rangeframe()+
  theme_tufte() +
  scale_x_continuous(breaks = seq(1, 11, by = 2)) +
  labs(x = "family size")

###################################################################################
################## there are only 3 genes HSPA13, PTPMT1 & PIK3C3 #################
reformed_para_cmb %>% 
  dplyr::filter(hgnc_symbol %in% c("HSPA13", "PTPMT1", "PIK3C3")) %>% 
  separate(para_name ,sep = "_", into = c("main_gene", "paralog_pairs")) -> main3_CYCLOPS

# for the paralog A
all_chr_gene %>% 
  dplyr::mutate(chr_pos = paste(chromosome_name, str_sub(band, 1,1), sep = "")) %>% 
  dplyr::select(hgnc_symbol, chr_pos)-> gene_chrpos

# for the paralog B
head(gene_chrpos)
gene_chrpos %>% 
  dplyr::rename(paralog_pairs = hgnc_symbol) -> paralog_chrpos
#colnames(gene_chrpos)[1] = "paralog_pairs"
main3_CYCLOPS %>% 
  left_join(gene_chrpos, by = "hgnc_symbol" ) %>% 
  left_join(paralog_chrpos, by ="paralog_pairs" ) %>% 
  dplyr::mutate(match = ifelse(chr_pos.x == chr_pos.y, 1, 0)) %>% 
  dplyr::filter(mean_identical_score >= 20) %>% 
  dplyr::select(hgnc_symbol,paralog_pairs, chr_pos.x, chr_pos.y) %>% View()


################ non-sig CYCLOPS #############
####### gene on the chromosome loss arm ##########
all_chr_gene %>% 
  dplyr::filter(!(chromosome_name %in% c("Y", "X"))) %>% 
  dplyr::mutate(chr_band = paste(chromosome_name, str_sub(band,1,1), sep = "")) %>% 
  dplyr::filter(chr_band %in% c("1p", "3p", "4p", "4q", "5q", "6p", "6q", "8p", "9p", "9q",
                                "10p", "10q", "11p", "11q", "12p", "13q", "14q", "15q", "16p",
                                "16q", "17p", "18p", "18q", "19p", "19q", "20p", "21q", "22q")) %>% 
  dplyr::select(hgnc_symbol) %>% 
  unlist(use.names = F) -> loss_chr_gene_lists

# non_sig_cyclops
non_sig_cyclops = loss_chr_gene_lists[-which(loss_chr_gene_lists %in% CYCLOPS_family_size$df_CYCLOPS.dep_gene)]

#family_size_sum
non_sig_cyclops_df = data.frame(non_sig_cyclops)
colnames(non_sig_cyclops_df) = "hgnc_symbol"
non_sig_cyclops_df %>% 
  left_join(family_size_sum,  by = "hgnc_symbol") %>% 
  dplyr::mutate(family_size_new = ifelse(is.na(family_size), 0, family_size)) %>% 
  dplyr::filter(family_size_new <=20) -> non_sig_cyclops_gene_family_size

non_sig_cyclops_gene_family_size %>% 
  dplyr::mutate(family_size_new = family_size_new + 1) %>% 
  dplyr::select(hgnc_symbol, family_size_new) -> non_sig_cyclops_gene_family_size_mdf

# plot
ggplot() +
  geom_density(
    data = CYCLOPS_family_size,
    aes(x = paralog_num),
    color = "#EF8536", fill = "#EF8536", alpha = 0.5
  ) +
  geom_density(
    data = non_sig_cyclops_gene_family_size,
    aes(x = family_size_new),
    color = "#7852A9", fill = "#7852A9", alpha = 0.5
  ) +
 
  geom_rangeframe() +
  theme_tufte() +
  labs(x = "Family Size", y = "Frequency")

 

# compare together to check the family size between sig and non-sig CYCLOPS based on the famliy size 
colnames(CYCLOPS_family_size_mdf) = c("hgnc_symbol", "Family_size")
CYCLOPS_family_size_mdf$condition = "Significant_CYCLOPS"
colnames(non_sig_cyclops_gene_family_size_mdf) = c("hgnc_symbol", "Family_size")
non_sig_cyclops_gene_family_size_mdf$condition = "Non-Sig_CYCLOPS"

# combine 2 tables 
CYCLOPS_cmb = rbind(CYCLOPS_family_size_mdf, non_sig_cyclops_gene_family_size_mdf)


# 300, 420
CYCLOPS_cmb %>% 
  dplyr::mutate(condition = factor(condition, levels = c("Significant_CYCLOPS", 
                                                         "Non-Sig_CYCLOPS"))) %>% 
  ggplot(aes(x = condition, y = Family_size)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10, hjust = 1, angle = 45)
  )


# barplot
#df_plot <- CYCLOPS_cmb %>% 
#  group_by(condition) %>% 
#  mutate(total = n()) %>% 
#  ungroup() %>% 
#  group_by(condition, Family_size) %>% 
#  dplyr::summarise(size_count = n(),
#            total = first(total),
#            .groups = "drop") %>% 
#  mutate(freq = size_count / total) %>% 
#  dplyr::mutate(condition = factor(condition, levels = c("Significant_CYCLOPS", "Non-Sig_CYCLOPS"))) %>% 
#  tidyr::complete(
#    condition,
#    Family_size,
#    fill = list(freq = 0, size_count = 0)
#  )

  
plot_df <- CYCLOPS_cmb %>%
  dplyr::count(condition, Family_size, name = "size_count") %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c("Significant_CYCLOPS", "Non-Sig_CYCLOPS")
    )
  ) %>%
  tidyr::complete(
    condition,
    Family_size,
    fill = list(size_count = 0)
  ) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(
    total = sum(size_count),
    freq = size_count / total
  ) %>%
  dplyr::ungroup()
  
  
# Plot the figure
ggplot(plot_df, aes(x = factor(Family_size), y = freq)) +
  
  geom_histogram(aes( fill = condition, color = condition), 
                 stat = "identity",
                 position = position_dodge(width = 0.9),
                 alpha= 0.8,width = 0.6) +
  geom_rangeframe() + 
  theme_tufte() +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(breaks = as.character(seq(1, 21, by = 2)), drop = FALSE) +
   labs(x = "Family size", y = "Frequency")

# chisqr test
plot_df %>% 
  dplyr::filter(Family_size == 1) -> df_plot_family_size1
  
p = 0.001263
chisq.test(
  matrix(
    c(18, 10,
      3770, 7433),
    nrow = 2,
    byrow = TRUE
  )
) 



# statistic test
#1.912e-06
t.test(CYCLOPS_cmb[CYCLOPS_cmb$condition == "Significant_CYCLOPS", ]$Family_size, 
       CYCLOPS_cmb[CYCLOPS_cmb$condition == "Non-Sig_CYCLOPS", ]$Family_size)



###################  Dependency difference between 2 groups ####################
# dependency filter 
# To judge that paralog could create dependency, we integrate the Depmap data and expression data from CCLE

# cell line info
metasheet_cellline = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Model.csv")

# filter the cell line are aneuploid
aneuploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")

aneuploid_dep %>% 
  dplyr::filter(Aneuploidy.score >= 7) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> aneuploid_cellline

df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")


# Gene dependency 24Q2 data set
df_gene_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /CRISPRGeneEffect.csv")

# change the name of data
colnames(df_gene_dep)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_Dep = sub("\\..*","",colnames(df_gene_dep)[2:dim(df_gene_dep)[2]])

# replace gene to colnames
colnames(df_gene_dep)[2:dim(df_gene_dep)[2]] = genelist_Dep
head(df_gene_dep[,1:6])

df_gene_dep_flt = df_gene_dep[df_gene_dep$DepMap_ID %in% aneuploid_cellline, ]
Gene_dependency_score = function(x){
  # get the gene
  gene = x
  # get the dependency
  gene_dep = df_gene_dep_flt[,which(colnames(df_gene_dep_flt) == gene)]
  
  if(!is_empty(gene_dep)){
    # calculate the median of the dependency
    median_dep = median(gene_dep[!is.na(gene_dep)])
    
    
  }else{
    median_dep = NA
  }
  
  return(median_dep)
  
}

# get the dependency
#dependency
sig_CYCLOPS_dep = lapply(CYCLOPS_family_size$df_CYCLOPS.dep_gene, Gene_dependency_score)
sig_CYCLOPS_dep_cmb = do.call(rbind, sig_CYCLOPS_dep)
sig_CYCLOPS_dep_cmb = as.data.frame(sig_CYCLOPS_dep_cmb)
sig_CYCLOPS_dep_cmb$condition = "sig"
colnames(sig_CYCLOPS_dep_cmb)[1] = "count"

# non_sig
nonsig_singleton_dep = lapply(non_sig_cyclops, Gene_dependency_score)
nonsig_singleton_dep_cmb = do.call(rbind, nonsig_singleton_dep)
nonsig_singleton_dep_cmb = as.data.frame(nonsig_singleton_dep_cmb)
nonsig_singleton_dep_cmb$condition = "nonsig"
colnames(nonsig_singleton_dep_cmb)[1] = "count"


# combine the file
dep_cmb_df = rbind(sig_CYCLOPS_dep_cmb, nonsig_singleton_dep_cmb)


dep_cmb_df %>% 
  dplyr::mutate(condition = factor(condition,levels = c("sig", "nonsig"))) %>% 
  ggplot(aes(x = condition,y = count)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(x = "", y = "dependency score") +
  ylim(c(-2,1))

#
 #p-value = 9.78e-15
wilcox.test(dep_cmb_df[dep_cmb_df$condition == "sig",]$count, 
            dep_cmb_df[dep_cmb_df$condition == "nonsig",]$count)



########### look at CYCLOPS and non-CYCLOPS genes 


################# CYCLOPS dependency analysis - neutral vs. chr loss for each classification ####################
CYCLOPS_family_size_mdf %>% 
  left_join(all_chr_gene, by = "hgnc_symbol") %>% 
  dplyr::select(hgnc_symbol, chromosome_name, band) %>% 
  dplyr::mutate(chr_name = paste("X",chromosome_name, str_sub(band ,1,1), sep  = "" )) -> sig_CYCLOPS_mdf

# aneuploidy cell aneuploid score condition  (only focus on aneuploidy cancer cell)
df_aneu_score %>% 
  dplyr::filter(X %in% aneuploid_cellline) -> aneuploid_status_flt

# function of normalised dependency 

SIG_CYCLOPS_dep = function(x){
  # get the gene name and chromsome postion
  
  df_tmp = sig_CYCLOPS_mdf[x,]
  
  # gene  
  
  gene_name = df_tmp[,1]
  
  # chromosome position
  
  chr_pos = df_tmp[,4]
  
  # get the certainc chromosome info
  loss_chr_condition = aneuploid_status_flt[,c(1, which(colnames(aneuploid_status_flt) == chr_pos))]
  #look a the certain chrosmome loss cell
  chr_loss_cells = loss_chr_condition[loss_chr_condition[,2] == -1, ]$X
  
  #look a the certain chrosmome neutral cell
  chr_neutral_cells = loss_chr_condition[loss_chr_condition[,2] == 0, ]$X
  
  # get the dependnecy score
  gene_dependency_flt = df_gene_dep[,c(1, which(colnames(df_gene_dep) ==gene_name ))]
  
  # get the loss dependency 
  gene_dependency_flt %>% 
    dplyr::filter(DepMap_ID %in% chr_loss_cells) %>% 
    .[,2] %>% unlist(use.names = F) -> loss_dependency
  
  # get the neutral dependency 
  gene_dependency_flt %>% 
    dplyr::filter(DepMap_ID %in% chr_neutral_cells) %>% 
    .[,2] %>% unlist(use.names = F) -> neutral_dependency
  
  # normalised the dependency
  median_dep_loss = median(loss_dependency[!is.na(loss_dependency)])
  
  median_dep_neutral = median(neutral_dependency[!is.na(neutral_dependency)])
  
  normm_loss_dep = median_dep_loss- median_dep_neutral
  
  norm_neutral_dep = median_dep_neutral- median_dep_neutral
  
  df_return = data.frame(gene_name,normm_loss_dep, norm_neutral_dep)
  
  return(df_return)
  
}

sig_CYCLOPS_dep_norm = lapply(c(1:dim(sig_CYCLOPS_mdf)[1]), SIG_CYCLOPS_dep)  
sig_CYCLOPS_dep_norm_cmb = do.call(rbind, sig_CYCLOPS_dep_norm)


sig_CYCLOPS_dep_norm_cmb %>% 
  gather(normm_loss_dep:norm_neutral_dep, key = "condition", value = "value") %>% 
  dplyr::mutate(condition = factor(condition, levels = c("normm_loss_dep", "norm_neutral_dep"))) -> sig_CYCLOPS_dep_mdf

sig_CYCLOPS_dep_mdf %>% 
  ggplot(aes(x = condition,y = value)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.2)) +
  coord_cartesian(ylim = c(-0.3, 0.3))+
  labs(x = "", y = "dependency score") -> P1


# Statistic
# p= 6.986e-12
wilcox.test(sig_CYCLOPS_dep_mdf[sig_CYCLOPS_dep_mdf$condition == "normm_loss_dep", ]$value, 
            sig_CYCLOPS_dep_mdf[sig_CYCLOPS_dep_mdf$condition == "norm_neutral_dep", ]$value)
################### non sig CYCLOPS ###############
non_sig_cyclops_gene_family_size_mdf %>% 
  left_join(all_chr_gene, by = "hgnc_symbol") %>% 
  dplyr::select(hgnc_symbol, chromosome_name, band) %>% 
  dplyr::mutate(chr_name = paste("X",chromosome_name, str_sub(band ,1,1), sep  = "" )) -> nonsig_CYCLOPS_mdf

# aneuploidy cell aneuploid score condition  (only focus on aneuploidy cancer cell)
df_aneu_score %>% 
  dplyr::filter(X %in% aneuploid_cellline) -> aneuploid_status_flt

# function of normalised dependency 

NONSIG_CYCLOPS_dep = function(x){
  # get the gene name and chromsome postion
  
  df_tmp = nonsig_CYCLOPS_mdf[x,]
  
  # gene  
  
  gene_name = df_tmp[,1]
 # print(gene_name)
  
  # chromosome position
  
  chr_pos = df_tmp[,4]
  
  # get the certainc chromosome info
  loss_chr_condition = aneuploid_status_flt[,c(1, which(colnames(aneuploid_status_flt) == chr_pos))]
  #look a the certain chrosmome loss cell
  chr_loss_cells = loss_chr_condition[loss_chr_condition[,2] == -1, ]$X
  
  #look a the certain chrosmome neutral cell
  chr_neutral_cells = loss_chr_condition[loss_chr_condition[,2] == 0, ]$X
  
  # get the dependnecy score
  gene_dependency_flt = df_gene_dep[,c(1, which(colnames(df_gene_dep) ==gene_name ))]
  if(!is.null(dim(gene_dependency_flt))){
    # get the loss dependency 
    gene_dependency_flt %>% 
      dplyr::filter(DepMap_ID %in% chr_loss_cells) %>% 
      .[,2] %>% unlist(use.names = F) -> loss_dependency
    
    # get the neutral dependency 
    gene_dependency_flt %>% 
      dplyr::filter(DepMap_ID %in% chr_neutral_cells) %>% 
      .[,2] %>% unlist(use.names = F) -> neutral_dependency
    
    # normalised the dependency
    median_dep_loss = median(loss_dependency[!is.na(loss_dependency)])
    
    median_dep_neutral = median(neutral_dependency[!is.na(neutral_dependency)])
    
    normm_loss_dep = median_dep_loss- median_dep_neutral
    
    norm_neutral_dep = median_dep_neutral- median_dep_neutral
  }else{
    normm_loss_dep = NA
    
    norm_neutral_dep = NA
  }
 
  
  df_return = data.frame(gene_name,normm_loss_dep, norm_neutral_dep)
  
  return(df_return)
  
}

nonsig_CYCLOPS_dep_norm = lapply(c(1:dim(nonsig_CYCLOPS_mdf)[1]), NONSIG_CYCLOPS_dep)  
nonsig_CYCLOPS_dep_norm_cmb = do.call(rbind, nonsig_CYCLOPS_dep_norm)

nonsig_CYCLOPS_dep_norm_cmb %>% 
  gather(normm_loss_dep:norm_neutral_dep, key = "condition", value = "value") %>% 
  dplyr::mutate(condition = factor(condition, levels = c("normm_loss_dep", "norm_neutral_dep"))) -> non_sig_CYCLOPS_dep_mdf
non_sig_CYCLOPS_dep_mdf %>% 
  ggplot(aes(x = condition,y = value)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(breaks = seq(-0.3, 0.3, by = 0.2)) +
  coord_cartesian(ylim = c(-0.3, 0.3))+
  labs(x = "", y = "dependency score") -> P2

# Statistic
# p= 0.1307
wilcox.test(non_sig_CYCLOPS_dep_mdf[non_sig_CYCLOPS_dep_mdf$condition == "normm_loss_dep", ]$value, 
            non_sig_CYCLOPS_dep_mdf[non_sig_CYCLOPS_dep_mdf$condition == "norm_neutral_dep", ]$value)

ggarrange(P1,P2, ncol = 2)

############## Check for the expression ##################

# retrieve gene expression data
CCLE_exp = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

# change the name of CCL2_exp data
colnames(CCLE_exp)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_exp = sub("\\..*","",colnames(CCLE_exp)[2:dim(CCLE_exp)[2]])

# replace gene to colnames
colnames(CCLE_exp)[2:dim(CCLE_exp)[2]] = genelist_exp

###########################################################
### sig CYCLOPS gene expression in aneuploid cancer cell###

SIG_CYCLOPS_exp = function(x){
  # get the gene name and chromsome postion
  
  df_tmp = sig_CYCLOPS_mdf[x,]
  
  # gene  
  
  gene_name = df_tmp[,1]
  
  # chromosome position
  
  chr_pos = df_tmp[,4]
  
  # get the certainc chromosome info
  loss_chr_condition = aneuploid_status_flt[,c(1, which(colnames(aneuploid_status_flt) == chr_pos))]
  #look a the certain chrosmome loss cell
  chr_loss_cells = loss_chr_condition[loss_chr_condition[,2] == -1, ]$X
  
  #look a the certain chrosmome neutral cell
  chr_neutral_cells = loss_chr_condition[loss_chr_condition[,2] == 0, ]$X
  
  # get the expression level
  gene_exp_flt = CCLE_exp[,c(1, which(colnames(CCLE_exp) ==gene_name ))]
  
  # get the loss expression 
  gene_exp_flt %>% 
    dplyr::filter(DepMap_ID %in% chr_loss_cells) %>% 
    .[,2] %>% unlist(use.names = F) -> loss_exp
  
  # get the neutral expression 
  gene_exp_flt %>% 
    dplyr::filter(DepMap_ID %in% chr_neutral_cells) %>% 
    .[,2] %>% unlist(use.names = F) -> neutral_exp
  
  # normalised the expression
  median_exp_loss = median(loss_exp[!is.na(loss_exp)])
  
  median_exp_neutral = median(neutral_exp[!is.na(neutral_exp)])
  
  normm_loss_exp= median_exp_loss- median_exp_neutral
  
  norm_neutral_exp = median_exp_neutral - median_exp_neutral
  
  df_return = data.frame(gene_name,normm_loss_exp, norm_neutral_exp)
  
  return(df_return)
  
}

sig_CYCLOPS_exp_norm = lapply(c(1:dim(sig_CYCLOPS_mdf)[1]), SIG_CYCLOPS_exp)  
sig_CYCLOPS_exp_norm_cmb = do.call(rbind, sig_CYCLOPS_exp_norm)


sig_CYCLOPS_exp_norm_cmb %>% 
  gather(normm_loss_exp:norm_neutral_exp, key = "condition", value = "value") %>% 
  dplyr::mutate(condition = factor(condition, levels = c("normm_loss_exp", "norm_neutral_exp"))) -> sig_CYCLOPS_exp_mdf

sig_CYCLOPS_exp_mdf %>% 
  ggplot(aes(x = condition,y = value)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(breaks = seq(-1, 0.75, by = 0.25)) +
  coord_cartesian(ylim = c(-1, 0.75)) +
  labs(x = "", y = "Relative expression") -> P3



################### non sig #########################
NONSIG_CYCLOPS_exp = function(x){
  # get the gene name and chromsome postion
  
  df_tmp = nonsig_CYCLOPS_mdf[x,]
  
  # gene  
  
  gene_name = df_tmp[,1]
  # print(gene_name)
  
  # chromosome position
  
  chr_pos = df_tmp[,4]
  
  # get the certainc chromosome info
  loss_chr_condition = aneuploid_status_flt[,c(1, which(colnames(aneuploid_status_flt) == chr_pos))]
  #look a the certain chrosmome loss cell
  chr_loss_cells = loss_chr_condition[loss_chr_condition[,2] == -1, ]$X
  
  #look a the certain chrosmome neutral cell
  chr_neutral_cells = loss_chr_condition[loss_chr_condition[,2] == 0, ]$X
  
  # get the exp score
  gene_exp_flt = CCLE_exp[,c(1, which(colnames(CCLE_exp) ==gene_name ))]
  if(!is.null(dim(gene_exp_flt))){
    # get the loss exp 
    gene_exp_flt %>% 
      dplyr::filter(DepMap_ID %in% chr_loss_cells) %>% 
      .[,2] %>% unlist(use.names = F) -> loss_exp
    
    # get the neutral exp 
    gene_exp_flt %>% 
      dplyr::filter(DepMap_ID %in% chr_neutral_cells) %>% 
      .[,2] %>% unlist(use.names = F) -> neutral_exp
    
    # normalised the exp
    median_exp_loss = median(loss_exp[!is.na(loss_exp)])
    
    median_exp_neutral = median(neutral_exp[!is.na(neutral_exp)])
    
    normm_loss_exp = median_exp_loss- median_exp_neutral
    
    norm_neutral_exp = median_exp_neutral- median_exp_neutral
  }else{
    normm_loss_exp = NA
    
    norm_neutral_exp = NA
  }
  
  
  df_return = data.frame(gene_name,normm_loss_exp, norm_neutral_exp)
  
  return(df_return)
  
}

nonsig_CYCLOPS_exp_norm = lapply(c(1:dim(nonsig_CYCLOPS_mdf)[1]), NONSIG_CYCLOPS_exp)  
nonsig_CYCLOPS_exp_norm_cmb = do.call(rbind, nonsig_CYCLOPS_exp_norm)

nonsig_CYCLOPS_exp_norm_cmb %>% 
  gather(normm_loss_exp:norm_neutral_exp, key = "condition", value = "value") %>% 
  dplyr::mutate(condition = factor(condition, levels = c("normm_loss_exp", "norm_neutral_exp"))) -> non_sig_CYCLOPS_exp_mdf
non_sig_CYCLOPS_exp_mdf %>% 
  ggplot(aes(x = condition,y = value)) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.25,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  scale_y_continuous(breaks = seq(-1, 0.75, by = 0.25)) +
  coord_cartesian(ylim = c(-1, 0.75)) +
  labs(x = "", y = "Relative expresssion") -> P4

# Statistic
# p= 4.905e-08
wilcox.test(non_sig_CYCLOPS_exp_mdf[non_sig_CYCLOPS_exp_mdf$condition == "normm_loss_exp", ]$value, 
            sig_CYCLOPS_exp_mdf[sig_CYCLOPS_exp_mdf$condition == "normm_loss_exp", ]$value)


ggarrange(P3,P4)
