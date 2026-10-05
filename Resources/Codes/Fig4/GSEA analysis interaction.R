# project: GSEA analysis of significant interaction targets 
# date: Mar.24.2025
# author: Yi

# library
library(tidyverse)
library(ggthemes)
library(clusterProfiler)
library(ggpubr)

# data
# load the data
sig_df = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/dep_summary_aneuploidy_chr_loss.xlsx")

# retrieve the significant interation genes
sig_df %>% 
  dplyr::filter(status == "interaction_dep") %>% 
  dplyr::select(dep_gene) %>% 
  unlist(use.names = F) ->sig_int_genes

# GSEA analysis from our GSEA analysis 
# function analysis 

GO_analysis = clusterProfiler::enrichGO(sig_int_genes,
                                        'org.Hs.eg.db', ont="BP" ,
                                        keyType = "SYMBOL",pvalueCutoff = 0.05)

View(GO_analysis[1:20,])
GO_analysis_df = as.data.frame(GO_analysis)
writexl::write_xlsx(GO_analysis_df, "Interactor_GO.xlsx")

# add the order ID and easy to pick up the annotaion
GO_analysis_df %>% 
  dplyr::mutate(ID_order = c(1:dim(GO_analysis_df)[1])) -> GO_item_order


GO_item_order %>% 
  .[c(1,2,3,4,5,6,7,9,27,38),] %>% 
  ggplot(aes(x = reorder(Description, -qvalue),y = -log10( qvalue))) +
  geom_point(aes(size = Count ,color = qvalue)) +
  coord_flip() +
  geom_rangeframe() +
  theme_tufte() +
  scale_color_gradient2(high ="#083D77", mid ="#ebebd3" ,low = "#E2725B",midpoint = 0.01)



