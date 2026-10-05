# project: comparing the CYCLOPS genes from our study with Rameen
# Date: Mar.21.2025
# Author: Yi

# library
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(ggvenn)

# load the data
sig_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")

rameen_target = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/CYCLOPS/Rameen_sig_genelist.xlsx")

# sig CYCLOPS target from my analysis
# the chormosome loss events 
sig_df %>% 
  dplyr::filter(status == "CYCLOPS") %>% 
  dplyr::select(dep_gene,aneuploid_loss_chr) -> CYCLOPS_df

# check the chromosome loss events in Rameen data
rameen_target %>% 
  dplyr::mutate(chromosome = paste("chr",Chromosome,sep = "")) %>%  
  dplyr::mutate(chr_loss = paste(chromosome, gsub("[^a-zA-Z]", "", Band), sep = "_")) %>% 
  dplyr::filter(chr_loss %in% unique(CYCLOPS_df$aneuploid_loss_chr) ) -> rameen_flt

intersect(CYCLOPS_df$dep_gene,rameen_flt$Gene)

# generate a VennDiagram
# list for 2 genesets
x = list(
  Our_list = CYCLOPS_df$dep_gene,
  CYCLOPS_paper = rameen_flt$Gene
)

#plot
ggvenn(
  x, 
  fill_color = c("#283593", "#C5CAE9"),
  stroke_size = 0.5, set_name_size = 4,show_percentage = F,fill_alpha = 0.75
)


# GSEA analysis from our GSEA analysis 
# function analysis 
library(clusterProfiler)

GO_analysis = clusterProfiler::enrichGO(CYCLOPS_df$dep_gene,
                          'org.Hs.eg.db', ont="BP" ,
                          keyType = "SYMBOL",pvalueCutoff = 0.05)

View(GO_analysis[1:20,])
GO_analysis_df = as.data.frame(GO_analysis)
writexl::write_xlsx(GO_analysis_df, "GO_CYCLOPS.xlsx")
# plot the GSEA analysis
library(MetBrewer)
GO_analysis_df %>% 
  .[c(1,2,3,5,6,7,8,10,12,13),] %>% 
  ggplot(aes(x = reorder(Description, -qvalue),y = -log10( qvalue))) +
  geom_point(aes(size = Count ,color = qvalue)) +
  coord_flip() +
  geom_rangeframe() +
  theme_tufte() +
  scale_color_gradient2(high ="#083D77", mid ="#ebebd3" ,low = "#E2725B",midpoint = 0.018)



