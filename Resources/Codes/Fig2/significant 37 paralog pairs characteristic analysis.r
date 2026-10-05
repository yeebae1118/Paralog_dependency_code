# Project: compare with significant 37 paralogs and non-significants pairs
# Date: Oct.27.2025
# Author: Yi

# Library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# load hg38 genomes
library(biomaRt)
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
  filter(!is.na(hsapiens_paralog_perc_id)) %>% 
  filter(hsapiens_paralog_associated_gene_name != "") %>% 
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
reformed_para_cmb = read.csv("paralog_pair_both_side.csv")
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


# In this analysis we only focus on the cell line are aneuploidy 
# load aneuploid score dataset

# CCLE aneuploid status annotation
aneuaploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")

# filter the cell line are aneuploid
aneuaploid_dep %>% 
  dplyr::filter(Aneuploidy.score >= 7) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> aneuploid_cellline

# get the aneuploidy condition
df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")

# filter out the naueploidy cells
df_aneu_score %>% 
  dplyr::filter(X %in% aneuploid_cellline) -> df_aneu_score_flt

# dependency filter 
# To judge that paralog could create dependency, we integrate the Depmap data and expression data from CCLE
# Gene dependency 24Q2 data set
df_gene_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /CRISPRGeneEffect.csv")

# cell line info
metasheet_cellline = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Model.csv")

# change the name of data
colnames(df_gene_dep)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_Dep = sub("\\.\\..*","",colnames(df_gene_dep)[2:dim(df_gene_dep)[2]])

# replace gene to colnames
colnames(df_gene_dep)[2:dim(df_gene_dep)[2]] = genelist_Dep
head(df_gene_dep[,1:6])

## analysis here: correlation and dep and exp
# retrieve gene expression data from CCLE
CCLE_exp = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

# change the name of CCL2_exp data
colnames(CCLE_exp)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
#dim(df_gene_dep)
genelist_exp = sub("\\.\\..*","",colnames(CCLE_exp)[2:dim(CCLE_exp)[2]])

# replace gene to colnames
colnames(CCLE_exp)[2:dim(CCLE_exp)[2]] = genelist_exp

# get the paralog pair dependency 
paralog_sum_with_status_aneu = readRDS("paralog_sum_with_status_aneu.RDS")

# read the file all significant 37 significant paralogs
sig_dep_df = readxl::read_xlsx("dep_summary_aneuploidy_chr_loss.xlsx")
sig_dep_df %>% 
  dplyr::filter(status == "Paralog")-> sig_37_paralog_df

# get into paralog_pair_flt
# get the paralog_pairs
sig_37_paralog_df %>% 
  dplyr::mutate(para_cmb_1 = paste(dep_gene,paralog_gene, sep = "_")) %>% 
  dplyr::mutate(para_cmb_2 = paste(paralog_gene,dep_gene, sep = "_")) %>% 
  dplyr::select(para_cmb_1,para_cmb_2) -> sig_37_paralog_cmb

# get the query paralog_combination
sig_sig_para_list = c(sig_37_paralog_cmb$para_cmb_1, sig_37_paralog_cmb$para_cmb_2)

# check for the paralog family size and identical score 
paralog_pair_flt %>% 
  dplyr::mutate(condition = ifelse(para_name %in% sig_sig_para_list, "sig", "non_sig")) %>% 
  dplyr::mutate(condition = factor(condition, levels = c("sig", "non_sig"))) -> family_identity_score_df

colors = c("#EF8536","#3A76AF")
# identity score 
family_identity_score_df %>% 
  ggplot(aes(x = condition, y = mean_identical_score)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  #geom_jitter( aes(color = status),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte()+
  theme(legend.position = "None")+
  labs(x = "", y = "Identical score") +
  scale_color_manual(values = c(colors[1], "#7852A9"))+
  scale_fill_manual(values = c(colors[1], "#7852A9")) -> P_identity

# perform test
# p-value = 0.003676
t.test(family_identity_score_df[family_identity_score_df$condition == "sig", ]$mean_identical_score, 
       family_identity_score_df[family_identity_score_df$condition == "non_sig", ]$mean_identical_score)


# family_size
family_identity_score_df %>% 
  ggplot(aes(x = condition, y = family_size)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  #geom_jitter( aes(color = status),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte()+
  theme(legend.position = "None")+
  labs(x = "", y = "family_size ") +
  scale_color_manual(values = c(colors[1], "#7852A9"))+
  scale_fill_manual(values = c(colors[1], "#7852A9")) -> P_family
  
#p-value = 2.64e-06
t.test(family_identity_score_df[family_identity_score_df$condition == "sig", ]$family_size, 
       family_identity_score_df[family_identity_score_df$condition == "non_sig", ]$family_size)


# combine the figure
ggarrange(P_identity,P_family,
          ncol = 2, 
          nrow = 1)


# identify the function of significants paralogs
sig_para_pairs = unique(c(sig_37_paralog_df$dep_gene, sig_37_paralog_df$paralog_gene))
library(clusterProfiler)
library(org.Hs.eg.db)
attributes <- c( "hgnc_symbol", "entrezgene_id")

# Fetch the gene information
gene_info_siglist <- getBM(attributes = attributes,
                           filters = "hgnc_symbol",
                           value = sig_para_pairs,
                           mart = mart,
                           useCache = T)

# GO analysis
ego <- enrichGO(gene          = sig_para_pairs,
                OrgDb         = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

View(ego[1:500,])

ego = as.data.frame(ego)
ego %>% 
  .[c(2,4,10,11,13,14, 18),] %>% 
  ggplot(aes(x = reorder(Description, -qvalue), y = -log10(qvalue))) +
  geom_point(aes(size = Count, color 
                 = qvalue), alpha = 0.75) +
  geom_rangeframe() +
  theme_tufte()+
  coord_flip()+
  scale_color_gradient2(low = "#E82035" ,mid = "#E5F3FD",  high = "#5371A9", midpoint = 5e-4)  +
  labs(x = "Description")
writexl::write_xlsx(ego, "sig_37_paralog_GO_analysis.xlsx")

# plot the GO analysis plot


################################################
## gene expression difference in sig_paralogs ##
################################################
# get the copy number dataset
# take the cell line from CCLE cnv data
CCLE_gene_cnv = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /PortalOmicsCNGeneLog2.csv")

# cell and patient info
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# rename the gene name on CCLE_gene_CNV
# change the name of data
colnames(CCLE_gene_cnv)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
dim(CCLE_gene_cnv)
genelist_Dep = sub("\\.\\..*","",colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]])

# replace gene to colnames
colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]] = genelist_Dep
head(CCLE_gene_cnv[,1:6])

# get the Y loss cell line ids
results_chrY <- getBM(attributes = c("chromosome_name", "entrezgene_id", "hgnc_symbol"),
                 filters = "chromosome_name", values = "Y", mart = mart)

# alternative 
gene_info %>% 
  dplyr::filter(chromosome_name == "Y") ->results_chrY

Y_chr_genename = results_chrY$hgnc_symbol[results_chrY$hgnc_symbol != ""]

# get the table as Y chromosome cnv
Y_cnv = CCLE_gene_cnv[,c(1,which(colnames(CCLE_gene_cnv) %in% Y_chr_genename))]

# get the male cell line
cell_line_id %>% 
  dplyr::filter(Gender == "male") %>% 
  dplyr::select(depMapID) %>% 
  unlist(use.names = F) -> male_cellline

# get the male Y chromomr cnv df
Y_cnv %>% 
  dplyr::filter(DepMap_ID %in% male_cellline) -> Y_cnv_male

# Y chr loss sample
Y_loss_sample = Y_cnv_male[rowSums(is.na(Y_cnv_male))>25,]$DepMap_ID 
# Y chr normal sample
Y_normal = Y_cnv_male$DepMap_ID [!(Y_cnv_male$DepMap_ID %in% Y_loss_sample)]

# comparing the expression of 37 genes on corresponding chromosomes
# get the cell loc genes and its corrpsodning chromosome position
sig_37_paralog_df %>% 
  dplyr::select(aneuploid_loss_chr, paralog_gene) -> paralog_loss_chr_gene_df
# get the chromosome location and apply for the loop
chr_pos = paralog_loss_chr_gene_df$aneuploid_loss_chr
chr_pos_edt = c()

for (chr in chr_pos) {
  if(chr == "chrY"){
    chr_pos_edt = c(chr_pos_edt, chr)
  }else{
    chr_tmp = str_sub(chr, 4, nchar(chr))
    chr_edt = paste("X", str_split(chr_tmp, pattern = "_")[[1]][1],str_split(chr_tmp, pattern = "_")[[1]][2],sep = "" )
    chr_pos_edt = c(chr_pos_edt, chr_edt)
  }
}

# assign the new chromosome position 

paralog_loss_chr_gene_df$chr_pos = chr_pos_edt

# compare the relative expression level 
#Sig relative paralog gene expression 
sig_para_relative_exp = function(x){
  
  # input list 
  df_tpm = paralog_loss_chr_gene_df[x,]
  
  # input paralog gene 
  gene = df_tpm$paralog_gene
  
  # input chromsome pos
  chr_pos = df_tpm$chr_pos
  
  if(chr_pos == "chrY"){
    
    # get the expression data
    exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene))]
    
    # loss exp 
    loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_loss_sample, 2]
    
    # median loss exp
    median_loss_exp = median(loss_exp_list)
    
    # neutral exp
    neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_normal, 2]
    
    median_neutral_exp = median(neutral_exp_list)
    
    # normalised the expression
    relative_loss_exp = median_loss_exp -median_neutral_exp
    relative_neutral_exp = median_neutral_exp -median_neutral_exp
    
  }else{
    # get the corresponding aneuploidy condition 
    df_chr_info_df = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == chr_pos))]
    # chr loss aneuploid cells
    chr_loss_sample = df_chr_info_df[df_chr_info_df[,2] == -1, ]$X
    
    # chr_neutral aneuploid cells
    chr_neutral_sample = df_chr_info_df[df_chr_info_df[,2] == 0, ]$X
    
    # get the expression data
    exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene))]
    
    # loss exp 
    loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_loss_sample, 2]
    
    # median loss exp
    median_loss_exp = median(loss_exp_list)
    
    # neutral exp
    neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_neutral_sample, 2]
    
    median_neutral_exp = median(neutral_exp_list)
    
    # normalised the expression
    relative_loss_exp = median_loss_exp -median_neutral_exp
    relative_neutral_exp = median_neutral_exp -median_neutral_exp
    
  }
  
  df_return = data.frame(relative_loss_exp, relative_neutral_exp)
  colnames(df_return) = c("chr_loss", "chr_neutral")
  
  return(df_return)
} 

sog_37_paralog_rel_exp = lapply(c(1:dim(paralog_loss_chr_gene_df)[1]), sig_para_relative_exp)
sog_37_paralog_rel_exp_cmb = do.call(rbind,sog_37_paralog_rel_exp )
sog_37_paralog_rel_exp_cmb %>% 
  gather(chr_loss:chr_neutral, key = "condition", value = "value") -> sig_mdf_exp
sig_mdf_exp %>% 
  ggplot(aes(x = condition, y = value))+
  geom_boxplot(width = 0.25, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.15) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel. paralog gene expression (TPM)")+
  coord_cartesian(ylim = c(-1, 0.2)) -> P_sig_paralog_exp

# statistic analysis

#p = 2.73e-13
wilcox.test(sig_mdf_exp[sig_mdf_exp$condition == "chr_loss", ]$value, sig_mdf_exp[sig_mdf_exp$condition == "chr_neutral", ]$value)

# ggarrange
ggarrange(P_sig_paralog_exp, ncol = 2)

###############################################
##### significant paralog pair dependency #####
###############################################

# We will compare significant paralog pairs to non-significant pairs 
# the paralog pairs
# do throught the significant paralog pairs

# get the chromosome location and apply for the loop
chr_pos_para = sig_37_paralog_df$aneuploid_loss_chr
chr_pos_para_edt = c()

for (chr in chr_pos_para) {
  if(chr == "chrY"){
    chr_pos_para_edt = c(chr_pos_para_edt, chr)
  }else{
    chr_tmp = str_sub(chr, 4, nchar(chr))
    chr_edt = paste("X", str_split(chr_tmp, pattern = "_")[[1]][1],str_split(chr_tmp, pattern = "_")[[1]][2],sep = "" )
    chr_pos_para_edt = c(chr_pos_para_edt, chr_edt)
  }
}

# for the dependency chromosome position
chr_pos_dep = sig_37_paralog_df$paralog_chr
chr_pos_dep_edt = c()

for (chr in chr_pos_dep) {
  if(chr %in% c("chrY", "chrX_p","chrX_q")){
    chr_pos_dep_edt = c(chr_pos_dep_edt, chr)
  }else{
    chr_tmp = str_sub(chr, 4, nchar(chr))
    chr_edt = paste("X", str_split(chr_tmp, pattern = "_")[[1]][1],str_split(chr_tmp, pattern = "_")[[1]][2],sep = "" )
    chr_pos_dep_edt = c(chr_pos_dep_edt, chr_edt)
  }
}

#assign the chromsosome position
sig_37_paralog_df$chr_loss_position = chr_pos_para_edt
sig_37_paralog_df$chr_dep_position = chr_pos_dep_edt

sig_dep_norm = function(x){
  
  # get the input data
  input_df = sig_37_paralog_df[x,]
  
  # Dependent gene
  dep_gene = input_df$dep_gene
  
  # Loss chromosome
  loss_chr = input_df$chr_loss_position
  
  if(loss_chr == "chrY"){
    # get the expression data
    dep_flt_df = df_gene_dep[,c(1,which(colnames(df_gene_dep) == dep_gene))]
    
    # loss exp 
    loss_dep_list = dep_flt_df[dep_flt_df$DepMap_ID %in% Y_loss_sample, 2]
    
    # median loss exp
    median_loss_dep = median(loss_dep_list[!is.na(loss_dep_list)])
    
    # neutral exp
    neutral_dep_list = dep_flt_df[dep_flt_df$DepMap_ID %in% Y_normal, 2]
    
    median_neutral_dep = median(neutral_dep_list[!is.na(neutral_dep_list)])
    
    # normalised the expression
    relative_loss_dep = median_loss_dep -median_neutral_dep
    relative_neutral_dep = median_neutral_dep -median_neutral_dep
  }else{
    # get the corresponding aneuploidy condition 
    df_chr_info_df = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == loss_chr))]
    # chr loss aneuploid cells
    chr_loss_sample = df_chr_info_df[df_chr_info_df[,2] == -1, ]$X
    
    # chr_neutral aneuploid cells
    chr_neutral_sample = df_chr_info_df[df_chr_info_df[,2] == 0, ]$X
    
    # get the expression data
    dep_flt_df = df_gene_dep[,c(1,which(colnames(df_gene_dep) == dep_gene))]
    
    # loss dep 
    loss_dep_list = dep_flt_df[dep_flt_df$DepMap_ID %in% chr_loss_sample, 2]
    
    # median loss dep
    median_loss_dep = median(loss_dep_list[!is.na(loss_dep_list)])
    
    # neutral dep
    neutral_dep_list = dep_flt_df[dep_flt_df$DepMap_ID %in% chr_neutral_sample, 2]
    
    median_neutral_dep = median(neutral_dep_list[!is.na(neutral_dep_list)])
    
    # normalised the dep
    relative_loss_dep = median_loss_dep -median_neutral_dep
    relative_neutral_dep = median_neutral_dep -median_neutral_dep
    
  }
  
  df_return = data.frame(relative_loss_dep, relative_neutral_dep)
  colnames(df_return) = c("chr_loss", "chr_neutral")
  
  return(df_return)
}


# get the relative significant loss dep
sig_norm_dep = lapply(c(1:dim(sig_37_paralog_df)[1]), sig_dep_norm)
sig_norm_dep_cmb = do.call(rbind,sig_norm_dep)

sig_norm_dep_cmb %>% 
  gather(chr_loss:chr_neutral, key = "condition", value = "value") -> sig_mdf_dep
sig_mdf_dep %>% 
  ggplot(aes(x = condition, y = value))+
  geom_boxplot(width = 0.3, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  #geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.15) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel. paralog dependency score")+
  coord_cartesian(ylim = c(-1, 0.2)) -> P_sig_paralog_dep
  

#p_value
#p-value = 2.683e-15
t.test(sig_mdf_dep[sig_mdf_dep$condition == "chr_loss", ]$value, sig_mdf_dep[sig_mdf_dep$condition == "chr_neutral", ]$value)

# go for the non-sig paralog pairs 
family_identity_score_df %>% 
  dplyr::filter(condition == "non_sig") %>% 
  separate(para_name, into = c("para_gene_1", "para_gene_2"), sep = "_")-> non_sig_para_pairs

# get each gene chromosome location : dep gene and paralog gene

gene_position <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "band"),
                  filters = "hgnc_symbol", 
                  values =unique(c(non_sig_para_pairs$para_gene_1, non_sig_para_pairs$para_gene_2)), 
                  mart = mart)

# alternative 
gene_info %>% 
  dplyr::filter(hgnc_symbol %in% unique(c(non_sig_para_pairs$para_gene_1, non_sig_para_pairs$para_gene_2))) %>% 
  dplyr::select(c("hgnc_symbol", "chromosome_name", "band")) -> gene_position

# edit the genes and combine with 
gene_position %>% 
  dplyr::filter(chromosome_name %in% c("X", "Y", c(1:22))) %>% 
  dplyr::filter(!is.na(band)) %>% 
  dplyr::mutate(chr_position = ifelse(chromosome_name %in% c("X","Y"), 
                                      paste("chr", chromosome_name, sep = ""), 
                                      paste("X", chromosome_name,  str_sub(band, 1,1), sep = ""))) %>% 
  dplyr::select(hgnc_symbol, chr_position) -> gene_position_mdf

# get para_gene_1
gene_position_mdf %>% 
  dplyr::rename(para_gene_1 = hgnc_symbol) %>% 
  distinct(para_gene_1, .keep_all = T)-> para_gene_1_position_mdf

# get_para_gene_2
gene_position_mdf %>% 
  dplyr::rename(para_gene_2 = hgnc_symbol,
                chr_para2_position = chr_position) %>% 
  distinct(para_gene_2, .keep_all = T)-> para_gene_2_position_mdf


# combine the data

non_sig_para_pairs %>% 
  left_join(para_gene_1_position_mdf, by = "para_gene_1") %>% 
  left_join(para_gene_2_position_mdf, by = "para_gene_2") %>% 
  dplyr::filter(!is.na(chr_position) ) %>% 
  dplyr::filter(!is.na(chr_para2_position) ) -> non_sig_para_query_tbl


# write function to test non_significant paralog dependency
non_sig_para_dep = function(x){
  
  # input the data
  input_df = non_sig_para_query_tbl[x,]
  
  #print(x)
  # paralog_gene_pairs
  # paralog_1 and its chromosome position 
  paralog_gene1 = input_df$para_gene_1
  gene1_pos = input_df$chr_position
  
  #paralog_2 and its chromosome position 
  paralog_gene2 = input_df$para_gene_2
  gene2_pos = input_df$chr_para2_position
  
  # check gene2 loss of aneuploidy and check the dependency on gene 1
  if(gene2_pos != "chrX"){
    if(gene2_pos == "chrY"){
      # get the expression data
      dep_flt_df_1 = df_gene_dep[,c(1,which(colnames(df_gene_dep) == paralog_gene1))]
      
      if(!is.null(dim(dep_flt_df_1))){
        # loss exp 
        loss_dep_list_1 = dep_flt_df_1[dep_flt_df_1$DepMap_ID %in% Y_loss_sample, 2]
        
        # median loss exp
        median_loss_dep_1 = median(loss_dep_list_1[!is.na(loss_dep_list_1)])
        
        # neutral exp
        neutral_dep_list_1 = dep_flt_df_1[dep_flt_df_1$DepMap_ID %in% Y_normal, 2]
        
        median_neutral_dep_1 = median(neutral_dep_list_1[!is.na(neutral_dep_list_1)])
        
        # normalised the expression
        relative_loss_dep_1 = median_loss_dep_1 - median_neutral_dep_1
        relative_neutral_dep_1 = median_neutral_dep_1 -median_neutral_dep_1
        
      }else{
        # normalised the expression
        relative_loss_dep_1 = NA
        relative_neutral_dep_1 = NA
      }
     
    }else{
      
      # get the corresponding aneuploidy condition 
      df_chr_info_df_1 = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == gene2_pos))]
      # chr loss aneuploid cells
      chr_loss_sample_1 = df_chr_info_df_1[df_chr_info_df_1[,2] == -1, ]$X
      
      # chr_neutral aneuploid cells
      chr_neutral_sample_1 = df_chr_info_df_1[df_chr_info_df_1[,2] == 0, ]$X
      
      # get the expression data
      dep_flt_df_1 = df_gene_dep[,c(1,which(colnames(df_gene_dep) == paralog_gene1))]
      
      if(!is.null(dim(dep_flt_df_1))){
        
        # loss dep 
        loss_dep_list_1 = dep_flt_df_1[dep_flt_df_1$DepMap_ID %in% chr_loss_sample_1, 2]
        
        # median loss dep
        median_loss_dep_1 = median(loss_dep_list_1[!is.na(loss_dep_list_1)])
        
        # neutral dep
        neutral_dep_list_1 = dep_flt_df_1[dep_flt_df_1$DepMap_ID %in% chr_neutral_sample_1, 2]
        
        median_neutral_dep_1 = median(neutral_dep_list_1[!is.na(neutral_dep_list_1)])
        
        # normalised the dep
        relative_loss_dep_1 = median_loss_dep_1 -median_neutral_dep_1
        relative_neutral_dep_1 = median_neutral_dep_1 -median_neutral_dep_1
        
        
      }else{
        # normalised the expression
        relative_loss_dep_1 = NA
        relative_neutral_dep_1 = NA
      }
      
     
      
    }
  }else{
    relative_loss_dep_1 = NA
    relative_neutral_dep_1 = NA
  }
  
  
  
  
  # check gene1 loss of aneuploidy and check the dependency on gene 2
  
  # check gene2 loss of aneuploidy and check the dependency on gene 1
  if(gene1_pos != "chrX"){
    if(gene1_pos == "chrY"){
      # get the expression data
      dep_flt_df_2 = df_gene_dep[,c(1,which(colnames(df_gene_dep) == paralog_gene2))]
      
      if(!is.null(dim(dep_flt_df_2))){
        # loss exp 
        loss_dep_list_2 = dep_flt_df_2[dep_flt_df_2$DepMap_ID %in% Y_loss_sample, 2]
        
        # median loss exp
        median_loss_dep_2 = median(loss_dep_list_2[!is.na(loss_dep_list_2)])
        
        # neutral exp
        neutral_dep_list_2 = dep_flt_df_2[dep_flt_df_2$DepMap_ID %in% Y_normal, 2]
        
        median_neutral_dep_2 = median(neutral_dep_list_2[!is.na(neutral_dep_list_2)])
        
        # normalised the expression
        relative_loss_dep_2 = median_loss_dep_2 - median_neutral_dep_2
        relative_neutral_dep_2 = median_neutral_dep_2 -median_neutral_dep_2
        
      }else{
        # normalised the expression
        relative_loss_dep_2 = NA
        relative_neutral_dep_2 = NA
      }
      
    }else{
      
      # get the corresponding aneuploidy condition 
      df_chr_info_df_2 = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == gene1_pos))]
      if(!is.null(dim(df_chr_info_df_2))){
        
        # chr loss aneuploid cells
        chr_loss_sample_2 = df_chr_info_df_2[df_chr_info_df_2[,2] == -1, ]$X
        
        # chr_neutral aneuploid cells
        chr_neutral_sample_2 = df_chr_info_df_2[df_chr_info_df_2[,2] == 0, ]$X
        
        # get the expression data
        dep_flt_df_2 = df_gene_dep[,c(1,which(colnames(df_gene_dep) == paralog_gene2))]
        
        if(!is.null(dim(dep_flt_df_2))){
          
          # loss dep 
          loss_dep_list_2 = dep_flt_df_2[dep_flt_df_2$DepMap_ID %in% chr_loss_sample_2, 2]
          
          # median loss dep
          median_loss_dep_2 = median(loss_dep_list_2[!is.na(loss_dep_list_2)])
          
          # neutral dep
          neutral_dep_list_2 = dep_flt_df_2[dep_flt_df_2$DepMap_ID %in% chr_neutral_sample_2, 2]
          
          median_neutral_dep_2 = median(neutral_dep_list_2[!is.na(neutral_dep_list_2)])
          
          # normalised the dep
          relative_loss_dep_2 = median_loss_dep_2 -median_neutral_dep_2
          relative_neutral_dep_2 = median_neutral_dep_2 -median_neutral_dep_2
          
          
        }else{
          # normalised the expression
          
          relative_loss_dep_2 = NA
          relative_neutral_dep_2 = NA
        }
        
        
      }else{
        # normalised the expression
        
        relative_loss_dep_2 = NA
        relative_neutral_dep_2 = NA
      }
      
      
    }
  }else{
    # normalised the expression
    relative_loss_dep_2 = NA
    relative_neutral_dep_2 = NA
  }
  
  return_1 = data.frame(relative_loss_dep_1, 
                        relative_neutral_dep_1)
  colnames(return_1) = c("chr_loss", "chr_neutral")
  retunr_2 = data.frame(relative_loss_dep_2, 
                          relative_neutral_dep_2)
  colnames(retunr_2) = c("chr_loss", "chr_neutral")
  return_df = rbind(return_1, retunr_2)
  
  return(return_df)
  
  
}


nonsig_pair_dep = lapply(c(1:dim(non_sig_para_query_tbl)[1]), non_sig_para_dep)

nonsig_pair_dep_cmb = do.call(rbind, nonsig_pair_dep)

View(nonsig_pair_dep_cmb)

nonsig_pair_dep_cmb %>% 
  gather(chr_loss:chr_neutral, key = "condition", value = "value") -> nonsig_pair_deps
nonsig_pair_deps %>% 
  ggplot(aes(x = condition, y = value)) +
  geom_boxplot(width = 0.3, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  #geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.15) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel. paralog dependency score")+
  coord_cartesian(ylim = c(-1, 0.2)) -> P_nonsig_paralog_dep

# plot sig and non-sig paraplog pairs
ggarrange(P_sig_paralog_dep,P_nonsig_paralog_dep, ncol = 2)

p_values <- c(
  t.test(sig_mdf_dep[sig_mdf_dep$condition == "chr_loss", ]$value, 
              sig_mdf_dep[sig_mdf_dep$condition == "chr_neutral", ]$value)$p.value,
  
  t.test(nonsig_pair_deps[nonsig_pair_deps$condition == "chr_loss", ]$value, 
              nonsig_pair_deps[nonsig_pair_deps$condition == "chr_neutral", ]$value)$p.value
  
)
# adjust the P value
# p_value = 2.681942e-09 1.000000e+00
p_adjusted_bonf <- p.adjust(p_values, method = "bonferroni") # Bonferroni correction



############################################
##   37 significant paralogs correlation  ##
############################################

sig_para_dep_exp_cor = function(x){
  
  # input list 
  df_tpm = sig_37_paralog_df[x,]
  
  # input paralog gene_dep
  gene_dep = df_tpm$dep_gene
  
  # input paralog_loss gene
  gene_loss = df_tpm$paralog_gene
  
  # get the expresssion of gene_loss
  exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene_loss))]

  # get the dependency data
  dep_flt_list = df_gene_dep[,c(1, which((colnames(df_gene_dep) ==gene_dep )))]
  
  if(x == 31){
    exp_flt_df %>% 
      left_join(dep_flt_list, by = "DepMap_ID") %>% 
      dplyr::filter(.[[2]] > 0) %>% 
      dplyr::filter(DepMap_ID %in% aneuploid_cellline) -> df_flt_exp_dep
  }else{
    # combine the data
    exp_flt_df %>% 
      left_join(dep_flt_list, by = "DepMap_ID") %>% 
      dplyr::filter(.[[2]] > 1) %>% 
      dplyr::filter(DepMap_ID %in% aneuploid_cellline) -> df_flt_exp_dep
    
  }
 
  cor_score = cor.test(df_flt_exp_dep[,2],df_flt_exp_dep[,3])$estimate
  df_return =data.frame(gene_dep,gene_loss ,cor_score)
  return(df_return)
} 

# apply for all the sig paralogs 
cor_exp_dep_list = lapply(c(1:dim(sig_37_paralog_df)[1]), sig_para_dep_exp_cor)

cor_exp_dep_list_cmb = do.call(rbind,cor_exp_dep_list)

cor_exp_dep_list_cmb %>% 
  ggplot(aes(x = cor_score, y = ..scaled..)) +
  geom_density(fill = "#5271A8", color = "#5271A8", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray49") +
  #geom_rangeframe() +
  theme_tufte() +
  xlim(c(-1,1)) +
  labs(x = "Correlation", y = "Scaled density") 


#################################################################################################
#### in 37 sig: Are significant paralogs more essential by themselves than non-sig paralogs? ####
#################################################################################################

# sig_37_dep_genes 
sig_dep_genes = unique(sig_37_paralog_df$dep_gene)

# get the non_sig_dep_genes

non_sig_dep_genes = unique(c(non_sig_para_query_tbl$para_gene_1, non_sig_para_query_tbl$para_gene_2))
# filter the genes in sig_37 genes
non_sig_dep_genes_flt = non_sig_dep_genes[!(non_sig_dep_genes%in% sig_dep_genes)]

# generate vector for the analysis 
gene_lists_dep = c()
dependency_score_median = c()
for (gene_name in c(sig_dep_genes, non_sig_dep_genes_flt)) {
  
  #print(gene_name)
  # get the dependency data 
  # only on aneuploidy cells
  dep_flt_list_flt = df_gene_dep[,c(1, which((colnames(df_gene_dep) == gene_name)))]
  if(!is.null(dim(dep_flt_list_flt))){
    dep_flt_list_flt %>% 
      dplyr::filter(DepMap_ID %in% aneuploid_cellline) %>% 
      .[,2] %>%  
      unlist(use.names = F) -> dep_score_flt
    
    median_dep_score = median(dep_score_flt[!is.na(dep_score_flt)])
  }else{
    median_dep_score = NA
  }
  
 
  # add to the list
  gene_lists_dep = c(gene_lists_dep, gene_name)
  dependency_score_median = c(dependency_score_median, median_dep_score)
  
  
}

# combine the dependency of sig and non-sig paralog genes

cmb_paralog_dep = data.frame(gene_lists_dep, dependency_score_median)
cmb_paralog_dep %>% 
  dplyr::mutate(condition = ifelse(gene_lists_dep %in%sig_dep_genes, "sig", "non_sig" )) -> essential_depenency_paralog

# plot the figures
essential_depenency_paralog %>% 
  dplyr::mutate(condition = factor(condition, levels = c("sig", "non_sig" ))) %>% 
  ggplot(aes(x = condition, y = dependency_score_median)) +
  geom_boxplot(width = 0.3, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  #geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.15) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel. paralog dependency score")+
  coord_cartesian(ylim = c(-2, 0.5))  -> rel_para_dependecy_comp

# test
# p-value < 2.2e-16
wilcox.test(essential_depenency_paralog[essential_depenency_paralog$condition == "sig", ]$dependency_score_median, 
            essential_depenency_paralog[essential_depenency_paralog$condition == "non_sig", ]$dependency_score_median)


ggarrange(rel_para_dependecy_comp, ncol = 2)



############################################################################################################
#### in 37 sig: #Do significant paralogs share more protein-protein interactions than non-sig paralogs? ####
############################################################################################################

library(rbioapi)
interactor_num = lapply(c(sig_dep_genes, non_sig_dep_genes_flt), function(x){
  # gene name
  gene_name = x
  
  proteins <- gene_name
  proteins_mapped_1 <- rba_string_map_ids(ids = proteins,
                                          species = 9606)
  
  proteins_mapped_query = proteins_mapped_1$stringId
  if(!is.null(proteins_mapped_query)){
    int_partners <- rba_string_interaction_partners(ids = proteins_mapped_query,
                                                    species = 9606,
                                                    required_score = 500,
                                                    limit = 100) 
    
    number_int = dim(int_partners[int_partners$score >= 0.9,])[1]
    
    return_df = data.frame(gene_name, number_int)
    
    return(return_df)
  }
  
  
})

interactor_num_cmb = do.call(rbind, interactor_num)
saveRDS(interactor_num_cmb, "interactor_num_cmb.rds")
interactor_num_cmb = readRDS("interactor_num_cmb.rds")
interactor_num_cmb %>% 
  dplyr::mutate(condition = ifelse(gene_name %in% sig_dep_genes, "sig", "non_sig" )) %>% 
  mutate(condition = factor(condition, levels = c("sig", "non_sig"))) -> int_num_para
int_num_para %>% 
  ggplot(aes(condition, number_int)) +
  
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.3,
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
  ) -> rel_para_int_numb


# test
# p-value 4.346e-15
wilcox.test(int_num_para[int_num_para$condition == "sig", ]$number_int, 
            int_num_para[int_num_para$condition == "non_sig", ]$number_int)


ggarrange(rel_para_int_numb, ncol = 2)





########## question for the significant paralog pairs has more shared partner than non-significant one ############
# get the function to check for the shared protein and protein interaction partner in shared paralog partners 
CO_SHARED_INTERACTOR = function(paralog_pair){
  # get the paralog pair list
  
  co_shared_percentage_list = c()
  for (pair in paralog_pair) {
   # get the gene name
    gene_list = str_split(pair, pattern = "_")[[1]]
     
    # pair1
    gene_1 = gene_list[1]
    
    # pair2
    gene_2 = gene_list[2]
    
    # get into the interation genes
    
    # gene_1 interaction
    
    proteins_mapped_1 <- rba_string_map_ids(ids = gene_1,
                                            species = 9606)
    
    proteins_mapped_query = proteins_mapped_1$stringId
    if(!is.null(proteins_mapped_query)){
      int_partners_1 <- rba_string_interaction_partners(ids = proteins_mapped_query,
                                                      species = 9606,
                                                      required_score = 500,
                                                      limit = 100) 
      
      int_partners_1_sig = int_partners_1[int_partners_1$score >= 0.7,]
      
      gene_1_interaction = int_partners_1_sig$preferredName_B
      
    }
    
    
    # gene 2 interaction 
    proteins_mapped_2 <- rba_string_map_ids(ids = gene_2,
                                            species = 9606)
    
    proteins_mapped_query_2 = proteins_mapped_2$stringId
    if(!is.null(proteins_mapped_query_2)){
      int_partners_2 <- rba_string_interaction_partners(ids = proteins_mapped_query_2,
                                                        species = 9606,
                                                        required_score = 500,
                                                        limit = 100) 
      
      int_partners_2_sig = int_partners_2[int_partners_2$score >= 0.7,]
      
      gene_2_interaction = int_partners_2_sig$preferredName_B
      
    }
    
    # get the interaction gene
    intersection =intersect(gene_1_interaction, gene_2_interaction)
    
    # get the percentage
    
    shared_percentage = (length(intersection) / length(gene_1_interaction) + length(intersection) / length(gene_2_interaction)) / 2
    
    # add into the list
    co_shared_percentage_list = c(co_shared_percentage_list, shared_percentage)
    
  }
  
  return(co_shared_percentage_list)
}  
#### look the table of significant paralog pairs #####

sig_37_paralog_df %>% 
  dplyr::mutate(para_pair_gene = paste(dep_gene,paralog_gene, sep = "_")) -> sig_37_paralog_df_mdf

# get the info of co_shared interaction partner

# get th sig paralog shard interaction partner
sig_co_location = CO_SHARED_INTERACTOR(sig_37_paralog_df_mdf$para_pair_gene)
sig_co_location_df = data.frame(sig_co_location)
colnames(sig_co_location_df) = "value"
sig_co_location_df$condition = "sig"
# get the non-significant pairs
family_identity_score_df %>% 
  dplyr::filter(condition == "non_sig") -> non_sig_pairs

# get the non-sig interatcion partner
non_sig_co_location = CO_SHARED_INTERACTOR(non_sig_pairs$para_name)

# 

non_sig_co_location_df = data.frame(non_sig_co_location)
colnames(non_sig_co_location_df) = "value"
non_sig_co_location_df$condition = "non-sig"

cmb_df_co_interaction = rbind(sig_co_location_df, non_sig_co_location_df)

saveRDS(cmb_df_co_interaction, "sig_37_shared_interaction_rds")
library(ggthemes)
cmb_df_co_interaction %>% 
  dplyr::mutate(condition = factor(condition, 
                                   levels = c("sig", "non-sig"))) %>% 
  ggplot(aes(x = condition, y =value )) +
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.3,
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
  ) -> paired_interction


# test
# p-value 0.007333
wilcox.test(cmb_df_co_interaction[cmb_df_co_interaction$condition == "sig", ]$value, 
            cmb_df_co_interaction[cmb_df_co_interaction$condition == "non-sig", ]$value)


ggarrange(paired_interction, ncol = 2)

############################################################################################################
####       in 37 sig: # Are significant paralogs more highly-conserved than non-sig paralogs?           ####
############################################################################################################
library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

calc_gene_conservation <- function(
    genes,
    bigwig,
    region = c("exon","gene_body","promoter"),
    promoter_upstream = 2000,
    promoter_downstream = 200,
    txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    orgdb = org.Hs.eg.db
){
  region <- match.arg(region)
  
  # Map symbols -> ENTREZ
  map <- AnnotationDbi::select(orgdb, keys = genes, keytype = "SYMBOL", columns = "ENTREZID")
  map <- map[!is.na(map$ENTREZID), c("SYMBOL","ENTREZID")]
  map <- map[!duplicated(map$SYMBOL), ]
  rownames(map) <- map$ENTREZID
  if (!nrow(map)) stop("None of the provided symbols mapped to ENTREZ IDs.")
  
  # Build GRanges per gene
  if (region == "exon") {
    ranges_list <- GenomicFeatures::exonsBy(txdb, by = "gene")[map$ENTREZID]
  } else if (region == "gene_body") {
    genes_txdb <- GenomicFeatures::genes(txdb)
    ranges_list <- GenomicRanges::GRangesList(genes_txdb[map$ENTREZID])
  } else { # promoter
    tx_by_gene <- GenomicFeatures::transcriptsBy(txdb, "gene")[map$ENTREZID]
    pick <- lapply(tx_by_gene, function(gr) if (length(gr)==0) GRanges() else gr[which.max(width(gr))])
    pick <- GRangesList(pick)
    prom <- promoters(unlist(pick, use.names = FALSE),
                      upstream = promoter_upstream, downstream = promoter_downstream)
    split_idx <- rep(names(pick), lengths(pick))
    ranges_list <- split(prom, split_idx)
  }
  
  # Drop empties
  ok <- lengths(ranges_list) > 0
  ranges_list <- ranges_list[ok]
  if (length(ranges_list) == 0) stop("No genomic ranges found for the requested region type.")
  
  # Normalize seqlevel style to UCSC (chr1, chr2, ...)
  # Do this on the unlisted ranges, then split again.
  ul <- unlist(ranges_list, use.names = FALSE)
  seqlevelsStyle(ul) <- "UCSC"   # <-- SAFE: value is a character(1)
  ranges_list <- split(ul, rep(names(ranges_list), lengths(ranges_list)))
  
  # Helper: weighted stats
  wmean <- function(scores, widths) {
    if (!length(scores)) return(NA_real_)
    sum(scores * widths, na.rm = TRUE) / sum(widths[!is.na(scores)], na.rm = TRUE)
  }
  wmedian <- function(scores, widths) {
    if (!length(scores)) return(NA_real_)
    o <- order(scores)
    s <- scores[o]; w <- widths[o]
    cw <- cumsum(w)/sum(w)
    s[which.max(cw >= 0.5)]
  }
  
  # Prepare holders
  res_mean <- setNames(rep(NA_real_, length(ranges_list)), names(ranges_list))
  res_median <- res_mean
  res_nbases <- setNames(integer(length(ranges_list)), names(ranges_list))
  
  # Work chromosome-by-chromosome (based on the actual ranges)
  all_ranges <- unlist(ranges_list, use.names = FALSE)
  chrs <- unique(as.character(seqnames(all_ranges)))
  
  for (chr in chrs) {
    # pick all gene ranges on this chr
    chr_ranges <- lapply(ranges_list, function(gr) gr[as.character(seqnames(gr)) == chr])
    chr_ranges <- chr_ranges[lengths(chr_ranges) > 0]
    if (!length(chr_ranges)) next
    
    merged_chr <- reduce(unlist(GRangesList(chr_ranges), use.names = FALSE))
    bw_chr <- import(bigwig, which = merged_chr)
    if (length(bw_chr) == 0) next
    
    # compute per gene
    for (gene_id in names(chr_ranges)) {
      ol <- findOverlaps(chr_ranges[[gene_id]], bw_chr, ignore.strand = TRUE)
      if (!length(ol)) next
      q <- chr_ranges[[gene_id]][queryHits(ol)]
      s <- bw_chr[subjectHits(ol)]
      intersected <- pintersect(q, s, ignore.strand = TRUE)
      w <- width(intersected)
      scores <- mcols(s)$score
      res_mean[gene_id]   <- wmean(scores, w)
      res_median[gene_id] <- wmedian(scores, w)
      res_nbases[gene_id] <- sum(w)
    }
  }
  
  kept_dt <- data.frame(ENTREZID = names(ranges_list),
                        gene = map[names(ranges_list), "SYMBOL"], row.names = NULL)
  
  out <- data.frame(
    gene = kept_dt$gene,
    n_bases_scored = as.integer(res_nbases),
    mean_score = as.numeric(res_mean),
    median_score = as.numeric(res_median),
    stringsAsFactors = FALSE
  )
  out[!is.na(out$mean_score) & out$n_bases_scored > 0, ]
}


# Exon-level conservation using PhyloP
#non_sig_dep_genes_flt_2 = non_sig_dep_genes_flt[-which(non_sig_dep_genes_flt %in% "HMGB1P1")]
map_tbl <- AnnotationDbi::select(org.Hs.eg.db, 
                             keys = c(sig_dep_genes, 
                                      non_sig_dep_genes_flt),
                             keytype = "SYMBOL", 
                             columns = "ENTREZID")
data_entriz_id = names(GenomicFeatures::exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene"))

# filter unmapped the genes
map_tbl %>% 
  dplyr::filter(ENTREZID %in% data_entriz_id) -> conser_gene_lists_flt

cons_all_paralog_exon <- calc_gene_conservation(
  genes = conser_gene_lists_flt$SYMBOL,
  bigwig = "~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /hg38.phyloP100way.bw",
  region = "exon"
)
head(cons_all_paralog_exon)
cons_all_paralog_exon %>% 
  dplyr::mutate(condition = ifelse(gene %in% sig_dep_genes, "sig", "non_sig" )) %>% 
  mutate(condition = factor(condition, levels = c("sig", "non_sig"))) -> conserv_df_plot
conserv_df_plot %>% 
  ggplot(aes(condition, median_score)) +
  
  geom_boxplot(
    aes(color = condition, fill = condition),
    width = 0.3,
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
  labs(x = "", y = "PhyloP score")+
  ylim(c(-1,2.5))  -> PhyloP_plot

#2.542e-05
wilcox.test(conserv_df_plot[conserv_df_plot$condition == "sig", ]$median_score, 
            conserv_df_plot[conserv_df_plot$condition == "non_sig", ]$median_score)

# plot the figures 
ggarrange(PhyloP_plot, ncol = 2)


############################################################################################################################################################
####  in 37 sig: #Do significant paralogs on lost chromosomes exhibit less dosage compensation compared to non-significant paralogs on lost chromosomes?####
############################################################################################################################################################
sig_37_paralog_df %>% 
 
  mutate(
    # pull captures per row
    m = str_match(aneuploid_loss_chr, "^chr(\\d+|X|Y)_([pq])$"),
    repl = paste0("X", m[,2], m[,3]),
    # keep chrY as-is; otherwise use replacement
    aneuploid_loss_chr = if_else(aneuploid_loss_chr == "chrY", aneuploid_loss_chr, repl)
  ) %>% 
  dplyr::select(paralog_gene,aneuploid_loss_chr) -> sig_loc_gene_df

# non sig loc gene and pos
non_sig_para_query_tbl %>% 
  dplyr::select(para_gene_1,chr_position) -> non_sig_query_1
colnames(non_sig_query_1) = c("gene_names", "chr_pos")
non_sig_para_query_tbl %>% 
  dplyr::select(para_gene_2,chr_para2_position) -> non_sig_query_2
colnames(non_sig_query_2) = c("gene_names", "chr_pos")

non_sig_query_list = rbind(non_sig_query_1, non_sig_query_2)
# can not include chrX loss
non_sig_query_list %>% 
  distinct(gene_names, .keep_all = T) %>% 
  dplyr::filter(!(gene_names %in% sig_loc_gene_df$paralog_gene)) %>% 
  dplyr::filter(chr_pos != "chrX") %>% 
  dplyr::filter(gene_names %in% colnames(CCLE_exp))-> non_sig_query_list_flt

# CCLE cell line with each chromsome arm gain or loss info
#df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")
#get all the aneuploidy cells
df_aneu_score %>% 
  dplyr::filter(X %in% aneuploid_cellline) -> aneuploid_cell_df_flt
#chromosome pos
chr_pos_list = colnames(aneuploid_cell_df_flt)[-1]

# check for the dosage compensation 
dosage_comp = function(x){
  
  # shape of data 
  return_df = lapply(1:dim(x)[1], function(a){
    # gene_name 
    print(a)
  
    gene_name = unlist(x[a,1], use.names = F)
    #chrom_pos
    chr_pos =unlist(x[a,2], use.names = F)
    
    if(chr_pos %in% chr_pos_list){
      # get the corresponding cell lines
      df_flt_anno = aneuploid_cell_df_flt[,c(1,which(colnames(aneuploid_cell_df_flt) == chr_pos))]
      
      # get the loss cells
      loss_cell_name = df_flt_anno[df_flt_anno[,2] == -1,]$X 
      
      # get the neutral cells
      neutral_cell_name = df_flt_anno[df_flt_anno[,2] == 0,]$X 
      
      # get expression
      exp_flt = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene_name))]
      #loss exp median
      exp_loss = exp_flt[exp_flt[,1] %in% loss_cell_name,][,2]
      # get the median exp
      median_loss_exp = median(exp_loss[!is.na(exp_loss)])
      # neutral exp median
      
      exp_neutral = exp_flt[exp_flt[,1]%in% neutral_cell_name,][,2]
      
      median_neutral_exp = median(exp_neutral[!is.na(exp_neutral)])
      
      result_return = data.frame(gene_name, median_loss_exp, median_neutral_exp)
      
      return(result_return)
      
    }else if(chr_pos == "chrY"){
      # get expression
      exp_flt = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene_name))]
      #loss exp median
      exp_loss = exp_flt[exp_flt[,1]%in% Y_loss_sample,][,2]
      # get the median exp
      median_loss_exp = median(exp_loss[!is.na(exp_loss)])
      # neutral exp median
      
      exp_neutral = exp_flt[exp_flt[,1]%in% Y_normal,][,2]
      
      median_neutral_exp = median(exp_neutral[!is.na(exp_neutral)])
      
      
    }else{
      
      median_loss_exp = NA
      median_neutral_exp= NA
      result_return = data.frame(gene_name, median_loss_exp, median_neutral_exp)
      
      return(result_return)
    }
    
    
  })
  
  
  return_df_cmb = do.call(rbind, return_df)
  
  return(return_df_cmb)
  
  
}

# sig_gene_exp
sig_gene_exp_df = dosage_comp(sig_loc_gene_df)
# non_sig gene exp
non_sig_gene_exp_df = dosage_comp(non_sig_query_list_flt)

# expression combination
exp_cmb_loc_genes = rbind(sig_gene_exp_df, non_sig_gene_exp_df)

exp_cmb_loc_genes %>% 
  filter(!is.na(median_loss_exp), 
         !is.na(median_neutral_exp)) %>% 
  mutate(norm_exp_para = median_loss_exp - median_neutral_exp,
         norm_exp_control = median_neutral_exp- median_neutral_exp) %>% 
  mutate(condition = ifelse(gene_name %in% sig_gene_exp_df$gene_name, "sig","non_sig")) %>% 
  mutate(condition = factor(condition, levels = c("sig", "non_sig"))) %>% 
  gather(norm_exp_para : norm_exp_control, key= "status", value = "value")-> exp_df_normalised
exp_df_normalised %>% 
  ggplot(aes(x = status, y = value))+
  geom_boxplot(
    aes(color = status, fill = status),
    width = 0.3,
    outlier.shape = NA, outliers = F,
    alpha = 0.4
  ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c( "#7852A9",colors[1])) +
  scale_color_manual(values = c("#7852A9",colors[1] )) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  ylab("Normalised exp")+
  ylim(c(-1,1)) +
  facet_wrap(~condition)-> sig_paralog_compen_plot


# test
# p-value / (2.275e-06)
wilcox.test(exp_df_normalised[exp_df_normalised$condition == "sig" &exp_df_normalised$status == "norm_exp_para"  ,]$value, 
            exp_df_normalised[exp_df_normalised$condition == "non_sig" & exp_df_normalised$status == "norm_exp_para" ,]$value)
#sig median
sig_paralog_median <- median(
  exp_df_normalised$value[
    exp_df_normalised$condition == "sig" &
      exp_df_normalised$status == "norm_exp_para"
  ],
  na.rm = TRUE
)


#non-sig median
nonsig_paralog_median <- median(
  exp_df_normalised$value[
    exp_df_normalised$condition == "non_sig" &
      exp_df_normalised$status == "norm_exp_para"
  ],
  na.rm = TRUE
)
print(paste("sig_paralog median:", sig_paralog_median))
print(paste("nonsig_paralog median:", nonsig_paralog_median))
# plot
ggarrange(sig_paralog_compen_plot, ncol = 2)

#######################################################
##############corum protein complex####################
#######################################################
## they should be in the same complex, incase they can compensate the function each other ##
#library(BioPlex)
corum_df <- read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Aneuploidy patient datasets/corum_humanComplexes.txt")


sig_37_paralog_df %>% 
  
  mutate(
    # pull captures per row
    m = str_match(aneuploid_loss_chr, "^chr(\\d+|X|Y)_([pq])$"),
    repl = paste0("X", m[,2], m[,3]),
    # keep chrY as-is; otherwise use replacement
    aneuploid_loss_chr = if_else(aneuploid_loss_chr == "chrY", aneuploid_loss_chr, repl)
  ) %>% 
  dplyr::select(dep_gene,aneuploid_loss_chr) -> sig_loc_gene_df

# non sig loc gene and pos
non_sig_para_query_tbl %>% 
  dplyr::select(para_gene_1,chr_position) -> non_sig_query_1
colnames(non_sig_query_1) = c("gene_names", "chr_pos")
non_sig_para_query_tbl %>% 
  dplyr::select(para_gene_2,chr_para2_position) -> non_sig_query_2
colnames(non_sig_query_2) = c("gene_names", "chr_pos")

non_sig_query_list = rbind(non_sig_query_1, non_sig_query_2)
# can not include chrX loss
non_sig_query_list %>% 
  distinct(gene_names, .keep_all = T) %>% 
  dplyr::filter(!(gene_names %in% sig_loc_gene_df$dep_gene)) %>% 
  dplyr::filter(chr_pos != "chrX") %>% 
  dplyr::filter(gene_names %in% colnames(CCLE_exp))-> non_sig_query_list_flt


# to compare the corum number difference 
CORUM_COMPLEX_NUMBER = function(x){
  # get the gene 
  gene = x
  
  # query for the corum number
  complex_df = corum_df[grepl(gene, corum_df$subunits_gene_name), ]
  
  if(is.null(dim(complex_df)[1])){
    complex_numb = 0
  }else{
    complex_numb = dim(complex_df)[1]
  }
  
  
  return(complex_numb)
}

#sig_vulnerab_corum_num
sig_para_corum = lapply(sig_loc_gene_df$dep_gene, CORUM_COMPLEX_NUMBER)
sig_para_corum_cmb = do.call(rbind, sig_para_corum)
sig_para_corum_cmb = as.data.frame(sig_para_corum_cmb)
sig_para_corum_cmb$condition = "sig"
colnames(sig_para_corum_cmb)[1] = "count"

# non_sig
nonsig_para_corum = lapply(non_sig_query_list_flt$gene_names, CORUM_COMPLEX_NUMBER)
nonsig_para_corum_cmb = do.call(rbind, nonsig_para_corum)
nonsig_para_corum_cmb = as.data.frame(nonsig_para_corum_cmb)
nonsig_para_corum_cmb$condition = "nonsig"
colnames(nonsig_para_corum_cmb)[1] = "count"

# combine the file
corum_cmb_df = rbind(sig_para_corum_cmb, nonsig_para_corum_cmb)
corum_cmb_df %>% 
  dplyr::mutate(condition = factor(condition,levels = c("sig", "nonsig"))) %>% 
  ggplot(aes(x = condition,y = count)) +
  geom_boxplot(
  aes(color = condition, fill = condition),
  width = 0.3,
  outlier.shape = NA, outliers = F,
  alpha = 0.4
 ) +
  geom_rangeframe() +
  theme_tufte() +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_y_sqrt()+
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(x = "", y = "Protein complex count")

#
# p-value = 0.0006717
t.test(corum_cmb_df[corum_cmb_df$condition == "sig",]$count, 
            corum_cmb_df[corum_cmb_df$condition == "nonsig",]$count)



################# CORUM protein complex ############################
# check the expresson on protein interacts partner expression

# RNA level
# sig_37_dep_genes 
# sig gene query list
sig_gene_query = data.frame(sig_37_paralog_df$dep_gene, sig_37_paralog_df$chr_loss_position)
colnames(sig_gene_query) = c("gene", "loss_chr_pos")

# get the non_sig_dep_genes
non_sig_gene_query = data.frame(c(non_sig_para_query_tbl$para_gene_1, non_sig_para_query_tbl$para_gene_2), c(non_sig_para_query_tbl$chr_para2_position, non_sig_para_query_tbl$chr_position))
colnames(non_sig_gene_query) = c("gene", "loss_chr_pos")
non_sig_gene_query %>% 
  dplyr::filter(!(gene %in% sig_gene_query$gene)) %>% 
  dplyr::filter(loss_chr_pos  %in% c("chrY", colnames(df_aneu_score_flt)[2:dim(df_aneu_score_flt)[2]]))-> non_sig_gene_query_flt

# get the protein corum data
# get protein omics data from CCLE
# get the protein MS data
protein_ms_df = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /harmonized_MS_CCLE_Gygi.csv")

# annotation of protein
uniprot_anno = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /uniprot_hugo_entrez_id_mapping.csv")

# get the colnames of proitein 
MS_protein_colname = colnames(protein_ms_df)[2:dim(protein_ms_df)[2]]

# make a query position vector
position_list = c()

# loop fpr the query protein position 
for (name in MS_protein_colname) {
  pos = which(uniprot_anno$UniprotID == name)
  
  position_list = c(position_list, pos)
}

# change the column names to the gene name 
colname_symbol = uniprot_anno[position_list,]$Symbol

# assign the gene to the column 
colnames(protein_ms_df)[2:dim(protein_ms_df)[2]] = colname_symbol

# check the expression 
# in the condition of chr loss and chr neutral
complex_exp_sig = function(x){
  
  # get the gene name
  query_df = sig_gene_query[x, ]
  
  # get the gene name
  gene_query = query_df$gene
  print(gene_query)
  
  # located chromosome
  chr_loc = query_df$loss_chr_pos
  
  # add the token
  pattern <- paste0("(^|[;,:[:space:]])(", paste(gene_query, collapse = "|"), ")([;,:[:space:]]|$)")
  
  # check if the CYCLOPS involves in any complex
  complex = corum_df[grepl(pattern, corum_df$subunits_gene_name),]
  
  # empty vector
  relative_loss_exp_list = c()
  relative_neutral_exp_list = c()
  # get the complex gene names 
  if(!(dim(complex)[1] == 0)){
    
    complex_gene_list = c()
    for (i in c(1:dim(complex)[1])) {
      complex_tpm = complex[i,]
      gene_names = str_split(pattern = ";",complex_tpm$subunits_gene_name)[[1]]
      gene_names_flt = gene_names[-which(gene_names == gene_query)]
      complex_gene_list = c(complex_gene_list, gene_names_flt)
      
    }
    
    
    if(chr_loc == "chrY"){
      
      for (gene in complex_gene_list) {
        # get expression data
        exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene))]
        
        # loss exp 
        loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_loss_sample, 2]
        
        # median loss exp
        median_loss_exp = median(loss_exp_list)
        
        # neutral exp
        neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_normal, 2]
        
        median_neutral_exp = median(neutral_exp_list)
        
        # normalised the expression
        relative_loss_exp = median_loss_exp -median_neutral_exp
        relative_neutral_exp = median_neutral_exp -median_neutral_exp
        
        # add to the empty vetcor
        relative_loss_exp_list = c(relative_loss_exp_list, relative_loss_exp)
        relative_neutral_exp_list = c(relative_neutral_exp_list, relative_neutral_exp)
      }
    }else{
      #chr_loc <- gsub("chr([0-9]+)_([pq])", "X\\1\\2", chr_loc)
      
      for(gene in complex_gene_list){
        # get the corresponding aneuploidy condition 
        df_chr_info_df = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == chr_loc))]
        # chr loss aneuploid cells
        chr_loss_sample = df_chr_info_df[df_chr_info_df[,2] == -1, ]$X
        
        # chr_neutral aneuploid cells
        chr_neutral_sample = df_chr_info_df[df_chr_info_df[,2] == 0, ]$X
        if(gene %in% colnames(CCLE_exp)){
          # get the expression data
          exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene))]
          
          # loss exp 
          loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_loss_sample, 2]
          
          # median loss exp
          median_loss_exp = median(loss_exp_list)
          
          # neutral exp
          neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_neutral_sample, 2]
          
          median_neutral_exp = median(neutral_exp_list)
        }else{
          median_loss_exp = NA
          median_neutral_exp = NA
        }
        
        
        # normalised the expression
        relative_loss_exp = median_loss_exp -median_neutral_exp
        relative_neutral_exp = median_neutral_exp -median_neutral_exp
        
        
        # add to the empty vetcor
        relative_loss_exp_list = c(relative_loss_exp_list, relative_loss_exp)
        relative_neutral_exp_list = c(relative_neutral_exp_list, relative_neutral_exp)
      }
      
    }
    
    
    
    
  }else{
    relative_loss_exp_list = c(relative_loss_exp_list, NA)
    relative_neutral_exp_list = c(relative_neutral_exp_list, NA)
  }
  
  
  return_df = data.frame(relative_loss_exp_list, relative_neutral_exp_list)
  colnames(return_df) = c("chr_loss", "chr_neutral")
  
  return(return_df)
  
}

corum_complex_exp_para_sig = lapply(c(1: dim(sig_gene_query)[1]), complex_exp_sig)  

corum_complex_exp_para_sig_cmb = do.call(rbind, corum_complex_exp_para_sig)
corum_complex_exp_para_sig_cmb %>% 
  gather(chr_loss:chr_neutral, key = "condition", value = "value") -> corum_gene_exp_sig
corum_gene_exp_sig %>% 
  ggplot(aes(x = condition, y = value))+
  geom_boxplot(width = 0.25, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.05) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel.gene expression") +
  ylim(c(-0.7,0.55))-> CORUM_gene_expression_sig

# p-value < 2.2e-16
wilcox.test(corum_gene_exp_sig[corum_gene_exp_sig$condition == "chr_loss", ]$value, corum_gene_exp_sig[corum_gene_exp_sig$condition == "chr_neutral", ]$value)

# at non-sig
# check the expression 
# in the condition of chr loss and chr neutral
complex_exp_nonsig = function(x){
  
  # get the gene name
  query_df = non_sig_gene_query_flt[x, ]
  
  # get the gene name
  gene_query = query_df$gene
  print(gene_query)
  
  # located chromosome
  chr_loc = query_df$loss_chr_pos
  
  # add the token
  pattern <- paste0("(^|[;,:[:space:]])(", paste(gene_query, collapse = "|"), ")([;,:[:space:]]|$)")
  
  # check if the CYCLOPS involves in any complex
  complex = corum_df[grepl(pattern, corum_df$subunits_gene_name),]
  
  # empty vector
  relative_loss_exp_list = c()
  relative_neutral_exp_list = c()
  # get the complex gene names 
  if(!(dim(complex)[1] == 0)){
    
    complex_gene_list = c()
    for (i in c(1:dim(complex)[1])) {
      complex_tpm = complex[i,]
      gene_names = str_split(pattern = ";",complex_tpm$subunits_gene_name)[[1]]
      gene_names_flt = gene_names[-which(gene_names == gene_query)]
      complex_gene_list = c(complex_gene_list, gene_names_flt)
      
    }
    if(!is_empty(complex_gene_list)){
      if(chr_loc == "chrY"){
        
        for (gene in complex_gene_list) {
          # get expression data
          
          exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == "PRP4K"))]
          if(!is.null(dim(exp_flt_df))){
            # loss exp 
            loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_loss_sample, 2]
            
            # median loss exp
            median_loss_exp = median(loss_exp_list)
            
            # neutral exp
            neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% Y_normal, 2]
            
            median_neutral_exp = median(neutral_exp_list)
            
            # normalised the expression
            relative_loss_exp = median_loss_exp -median_neutral_exp
            relative_neutral_exp = median_neutral_exp -median_neutral_exp
          }else{
            relative_loss_exp = NA
            relative_neutral_exp = NA
          }
         
          
          # add to the empty vetcor
          relative_loss_exp_list = c(relative_loss_exp_list, relative_loss_exp)
          relative_neutral_exp_list = c(relative_neutral_exp_list, relative_neutral_exp)
        }
      }else{
        #chr_loc <- gsub("chr([0-9]+)_([pq])", "X\\1\\2", chr_loc)
        
        for(gene in complex_gene_list){
          # get the corresponding aneuploidy condition 
          df_chr_info_df = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == chr_loc))]
          # chr loss aneuploid cells
          chr_loss_sample = df_chr_info_df[df_chr_info_df[,2] == -1, ]$X
          
          # chr_neutral aneuploid cells
          chr_neutral_sample = df_chr_info_df[df_chr_info_df[,2] == 0, ]$X
          if(gene %in% colnames(CCLE_exp)){
            # get the expression data
            exp_flt_df = CCLE_exp[,c(1,which(colnames(CCLE_exp) == gene))]
            
            # loss exp 
            loss_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_loss_sample, 2]
            
            # median loss exp
            median_loss_exp = median(loss_exp_list)
            
            # neutral exp
            neutral_exp_list = exp_flt_df[exp_flt_df$DepMap_ID %in% chr_neutral_sample, 2]
            
            median_neutral_exp = median(neutral_exp_list)
          }else{
            median_loss_exp = NA
            median_neutral_exp = NA
          }
          
          
          # normalised the expression
          relative_loss_exp = median_loss_exp -median_neutral_exp
          relative_neutral_exp = median_neutral_exp -median_neutral_exp
          
          
          # add to the empty vetcor
          relative_loss_exp_list = c(relative_loss_exp_list, relative_loss_exp)
          relative_neutral_exp_list = c(relative_neutral_exp_list, relative_neutral_exp)
        }
        
      }
      
      
    }else{
      relative_loss_exp_list = c(relative_loss_exp_list, NA)
      relative_neutral_exp_list = c(relative_neutral_exp_list, NA)
    }
    
    
  }else{
    relative_loss_exp_list = c(relative_loss_exp_list, NA)
    relative_neutral_exp_list = c(relative_neutral_exp_list, NA)
  }
  
  
  return_df = data.frame(relative_loss_exp_list, relative_neutral_exp_list)
  colnames(return_df) = c("chr_loss", "chr_neutral")
  
  return(return_df)
  
}

corum_complex_exp_para_nonsig = lapply(c(1: dim(non_sig_gene_query_flt)[1]), complex_exp_nonsig)  

corum_complex_exp_para_nonsig_cmb = do.call(rbind, corum_complex_exp_para_nonsig)
corum_complex_exp_para_nonsig_cmb %>% 
  gather(chr_loss:chr_neutral, key = "condition", value = "value") -> corum_gene_exp_nonsig
corum_gene_exp_nonsig %>% 
  ggplot(aes(x = condition, y = value))+
  geom_boxplot(width = 0.25, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  #geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.05) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel.gene expression") +
  ylim(c(-0.7,1))-> CORUM_gene_expression_nonsig

wilcox.test(corum_gene_exp_nonsig[corum_gene_exp_nonsig$condition == "chr_loss", ]$value, corum_gene_exp_nonsig[corum_gene_exp_nonsig$condition == "chr_neutral", ]$value)


########################at the protein level#####################
complex_protein_sig = function(x){
  
  # get the gene name
  query_df = sig_gene_query[x, ]
  
  # get the gene name
  gene_query = query_df$gene
  print(gene_query)
  
  # located chromosome
  chr_loc = query_df$loss_chr_pos
  
  # add the token
  pattern <- paste0("(^|[;,:[:space:]])(", paste(gene_query, collapse = "|"), ")([;,:[:space:]]|$)")
  
  # check if the CYCLOPS involves in any complex
  complex = corum_df[grepl(pattern, corum_df$subunits_gene_name),]
  
  # empty vector
  relative_loss_exp_list = c()
  relative_neutral_exp_list = c()
  # get the complex gene names 
  if(!(dim(complex)[1] == 0)){
    
    complex_gene_list = c()
    for (i in c(1:dim(complex)[1])) {
      complex_tpm = complex[i,]
      gene_names = str_split(pattern = ";",complex_tpm$subunits_gene_name)[[1]]
      gene_names_flt = gene_names[-which(gene_names == gene_query)]
      complex_gene_list = c(complex_gene_list, gene_names_flt)
      
    }
    
    
    if(chr_loc == "chrY"){
      
      for (gene in complex_gene_list) {
        # get expression data
        
        exp_flt_df = protein_ms_df[,c(1,which(colnames(protein_ms_df) == gene))]
        if(!is.null(dim(exp_flt_df))){
          # loss exp 
          loss_exp_list = exp_flt_df[exp_flt_df$X %in% Y_loss_sample, 2]
          
          # median loss exp
          median_loss_exp = median(loss_exp_list)
          
          # neutral exp
          neutral_exp_list = exp_flt_df[exp_flt_df$X %in% Y_normal, 2]
          
          median_neutral_exp = median(neutral_exp_list)
          
          # normalised the expression
          relative_loss_exp = median_loss_exp -median_neutral_exp
          relative_neutral_exp = median_neutral_exp -median_neutral_exp
        }else{
          relative_loss_exp = NA
          relative_neutral_exp = NA
        }
        
        
        # add to the empty vetcor
        relative_loss_exp_list = c(relative_loss_exp_list, relative_loss_exp)
        relative_neutral_exp_list = c(relative_neutral_exp_list, relative_neutral_exp)
      }
    }else{
      #chr_loc <- gsub("chr([0-9]+)_([pq])", "X\\1\\2", chr_loc)
      
      for(gene in complex_gene_list){
        # get the corresponding aneuploidy condition 
        df_chr_info_df = df_aneu_score_flt[,c(1, which(colnames(df_aneu_score_flt) == chr_loc))]
        # chr loss aneuploid cells
        chr_loss_sample = df_chr_info_df[df_chr_info_df[,2] == -1, ]$X
        
        # chr_neutral aneuploid cells
        chr_neutral_sample = df_chr_info_df[df_chr_info_df[,2] == 0, ]$X
        if(gene %in% colnames(protein_ms_df)){
          # get the expression data
          exp_flt_df = protein_ms_df[,c(1,which(colnames(protein_ms_df) == gene))]
          
          # loss exp 
          loss_exp_list = exp_flt_df[exp_flt_df$X %in% chr_loss_sample, 2]
          
          # median loss exp
          median_loss_exp = median(loss_exp_list)
          
          # neutral exp
          neutral_exp_list = exp_flt_df[exp_flt_df$X %in% chr_neutral_sample, 2]
          
          median_neutral_exp = median(neutral_exp_list)
        }else{
          median_loss_exp = NA
          median_neutral_exp = NA
        }
        
        
        # normalised the expression
        relative_loss_exp = median_loss_exp -median_neutral_exp
        relative_neutral_exp = median_neutral_exp -median_neutral_exp
        
        
        # add to the empty vetcor
        relative_loss_exp_list = c(relative_loss_exp_list, relative_loss_exp)
        relative_neutral_exp_list = c(relative_neutral_exp_list, relative_neutral_exp)
      }
      
    }
    
    
    
    
  }else{
    relative_loss_exp_list = c(relative_loss_exp_list, NA)
    relative_neutral_exp_list = c(relative_neutral_exp_list, NA)
  }
  
  
  return_df = data.frame(relative_loss_exp_list, relative_neutral_exp_list)
  colnames(return_df) = c("chr_loss", "chr_neutral")
  
  return(return_df)
  
}

corum_complex_protein_para_sig = lapply(c(1: dim(sig_gene_query)[1]), complex_protein_sig)  

corum_complex_protein_para_sig_cmb = do.call(rbind, corum_complex_protein_para_sig)
corum_complex_protein_para_sig_cmb %>% 
  gather(chr_loss:chr_neutral, key = "condition", value = "value") -> corum_gene_protein_sig
corum_gene_protein_sig %>% 
  ggplot(aes(x = condition, y = value))+
  geom_boxplot(width = 0.25, 
               aes(color = condition, fill = condition),
               alpha = 0.5,
               outliers = F) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.05) +
  geom_rangeframe()+
  theme_tufte() +
  theme(legend.position = "none") +
  scale_color_manual(values = c(colors[1], "#7852A9")) +
  scale_fill_manual(values = c(colors[1], "#7852A9")) +
  labs(x = "", y = "Rel.gene expression") +
  ylim(c(-0.7,0.55))-> CORUM_gene_protein_sig

wilcox.test(corum_gene_protein_sig[corum_gene_protein_sig$condition == "chr_loss", ]$value, corum_gene_protein_sig[corum_gene_protein_sig$condition == "chr_neutral", ]$value)


# 350 400
ggarrange(CORUM_gene_expression_sig, CORUM_gene_protein_sig)
