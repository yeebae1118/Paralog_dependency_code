# Project: 72 significant interactors and correlation
# Date: Nov.03.2025
# Author: Yi 


# library
library(tidyverse)
library(ggthemes)
library(ggpubr)
#library(biomaRt)

# load data
# dataset
# CCLE aneuploid status annotation
aneuploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")
# CCLE cell line with each chromsome arm gain or loss info
df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")
#cell line annotation
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# get aneuploid cells
# filter the cell line are aneuploid
aneuploid_dep %>% 
  dplyr::filter(Aneuploidy.score >= 7) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> aneuploid_cellline

# get the significant partner
# read the file all significant 37 significant paralogs
sig_dep_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")


# get the cell interactors dependency
sig_dep_df %>% 
  dplyr::filter(status == "interaction_dep") -> int_df_flt

# load hg38 genomes
library(biomaRt)
mart <- useMart(biomart="ensembl",
                dataset="hsapiens_gene_ensembl")
attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position",
                "end_position", "strand", "gene_biotype", "description","band")

# Fetch the gene information
gene_info <- getBM(attributes = attributes, mart = mart,useCache = TRUE)
# filter the chromosome 
gene_info %>% 
  dplyr::filter(chromosome_name %in% c(1:22,"Y","X")) %>% 
  dplyr::filter(hgnc_symbol != "") %>% 
  distinct(hgnc_symbol, .keep_all = T) %>% 
  dplyr::filter(gene_biotype == "protein_coding") %>% 
  dplyr::mutate(chr_pos = ifelse(chromosome_name %in% c("X", "Y"), paste("chr", chromosome_name, sep = ""),
                                 paste(paste("chr",chromosome_name, sep = "" ), str_sub(band, 1,1),sep = "_"))) %>% 
  dplyr::select(hgnc_symbol, chr_pos)-> all_chr_gene
# library
library(rbioapi)

# write function to get the gene interactors
interactor_number_likelihood = function(x){
  # input the data 
  input_df =  int_df_flt[x,]
  
  # dep_gene
  dep_gene = input_df$dep_gene
  
  # loss chromosome
  loss_chr_spec = input_df$aneuploid_loss_chr
  
  # get the dep gene 
  proteins_mapped_1 <- rba_string_map_ids(ids = dep_gene,
                                          species = 9606)
  
  
  proteins_mapped_query = proteins_mapped_1$stringId
  
  int_partners <- rba_string_interaction_partners(ids = proteins_mapped_query,
                                                  species = 9606,
                                                  required_score = 500,
                                                 limit = 200) # set the limit interactor number 200
  # cutoff at 0.7 and assign the chromosome position
  int_partners %>% 
    dplyr::filter(score > 0.5) %>% 
    dplyr::rename( hgnc_symbol = preferredName_B) %>% 
    left_join(all_chr_gene, by = "hgnc_symbol") -> df_int_partner
  
  # get the summary of number of interactors
  interactor_count = dim(df_int_partner)[1]
  df_int_partner %>% 
    dplyr::filter(chr_pos == loss_chr_spec) -> inter_loss_chr_flt
  
  # get the count the interactors 
  if(!is.null(dim(inter_loss_chr_flt))){
    interactor_loss_count = dim(inter_loss_chr_flt)[1]
  }else{
    interactor_loss_count = 0
  }
  
  
  return_df = data.frame(dep_gene, interactor_count, interactor_loss_count)
  colnames(return_df) = c("gene", "total_inter_count", "loss_chr_inter_count")
  
  return(return_df)
 
    
}

# lapply the function to the data
summary_inter_73 = lapply(c(1:dim(int_df_flt)[1]), interactor_number_likelihood)

summary_inter_73_cmb = do.call(rbind, summary_inter_73)

summary_inter_73_cmb %>% 
  dplyr::mutate(inter_count = loss_chr_inter_count) %>% 
  group_by(inter_count) %>% 
  dplyr::mutate(total_inter = sum(total_inter_count)) %>% 
  dplyr::mutate(specific_inter = sum(loss_chr_inter_count)) %>% 
  ungroup() %>% 
  distinct(inter_count, .keep_all = T) %>% 
  dplyr::mutate(freq_inter = specific_inter/total_inter )  -> df_tpm

#weigth
library(MetBrewer)
df_tpm$pc = predict(prcomp(~inter_count+freq_inter, df_tpm))[,1]
df_tpm %>% 
  ggplot(aes(x = inter_count,y = freq_inter)) +
  geom_point(aes(color = pc,,size = freq_inter),alpha = 0.8) +
  geom_smooth(method = "lm",alpha = 0.1) +
  geom_rangeframe() +
  theme_tufte() + 
  scale_color_gradient2(low = "#FCAE1E" ,mid = "#E2DDFC",  high = "#543D7B") + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(1, 16, by = 4)) +
  labs(x = "number of interator partner on one chromosome", y = "possibility of significant dependent targets")

#8.476e-06
#  cor 
#0.9345387 
cor_result = cor.test(df_tpm$inter_count, df_tpm$freq_inter)
  
# WE LOOK AT HYPK 
  
# get the dep gene 
proteins_mapped_HYPK <- rba_string_map_ids(ids = "HYPK",
                                        species = 9606)


HYPK_mapped_query = proteins_mapped_HYPK$stringId

HYPK_int_partners <- rba_string_interaction_partners(ids = HYPK_mapped_query,
                                                species = 9606,
                                                required_score = 500,
                                                limit = 100) 
HYPK_int_partners$preferredName_B

# library for ggraph
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)
library(scatterpie)

# read HYPK interaction network
HYPK_network = read.delim("HYPK_string_interactions.tsv")
colnames(all_chr_gene)[1] = "X.node1"

HYPK_network %>% 
  left_join(all_chr_gene, by = "X.node1") %>% 
  dplyr::select(X.node1, node2, chr_pos) -> HYPK_network_edt

unique(HYPK_network_edt$chr_pos)
# 
g <- graph_from_data_frame(
  d = HYPK_network_edt[,1:2],      # edges
  directed = FALSE       # or TRUE if directional
)

g <- igraph::simplify(g)

HYPK_network_edt %>% 
  distinct(X.node1, .keep_all = T) -> network_dist

V(g)$grp <- network_dist$chr_pos
V(g)$grp <- factor(V(g)$grp, levels = unique(network_dist$chr_pos))
library(oaqc)
bb <- layout_as_backbone(g, keep = 0.4)
E(g)$col <- F
E(g)$col[bb$backbone] <- T

#
library(MetBrewer)

nodes_to_label <- c("HYPK", "NAA15", "NAA11")

# plot
ggraph(g,
       layout = "manual",
       x = bb$xy[, 1],
       y = bb$xy[, 2]) +
  geom_edge_link0(aes(col = col), width = 0.2) +
  geom_node_point(aes(fill = grp), shape = 21, size = 3) +
  geom_node_text(
    aes(label = ifelse(name %in% nodes_to_label, name, "")),
    repel = TRUE,
    size = 3
  ) +
  geom_mark_hull(
    aes(x, y, group = grp, fill = grp, label=grp),
    concavity = 4,
    expand = unit(2, "mm"),
    alpha = 0.25
  ) +
  
 
  scale_color_manual(values  = c(rep("#0000004D", time = 2), "#e58068","#0000004D" , "#6b7cb9", rep("#0000004D", time = 10))) +
  scale_fill_manual(values  = c(rep("#0000004D", time = 2), "#e58068", "#0000004D", "#6b7cb9", rep("#0000004D", time = 10))) +
  scale_edge_color_manual(values = c(rgb(0, 0, 0, 0.3), rgb(0, 0, 0, 1))) +
  theme_graph()+
  theme(legend.position = "none")
