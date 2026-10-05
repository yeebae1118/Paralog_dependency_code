# checking the aneuploidy score in 8p loss cells
# date : Oct 28.2024
# Author: Yi

# library
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(MetBrewer)

# dataset
# CCLE aneuploid status annotation
aneuploid_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/aneuploidy_scores.csv")
# CCLE cell line with each chromsome arm gain or loss info
df_aneu_score = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/arm_call_scores.csv")
#cell line annotation
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# color set
colors = c("#EF8536","#3A76AF")


# definition: 
# Aneuploidy score  median = 25: high aneuploidy
# Near-euploid; median=3
# to avoid the bias from other chromosome and get accuracy dependecy on 7p gain
# we have to set up an aneuploidy score cut-off

#get the 8p loss cells
df_aneu_score %>% 
  filter(X8p == -1) %>% 
  dplyr::select(X) %>%
  unlist(use.names = F) -> chr_8p_loss

# median aneuploidy for chr8p loss
aneuploid_dep %>% 
  dplyr::filter(DepMap_ID  %in% chr_8p_loss) -> chr8p_loss_flt
median(chr8p_loss_flt$Aneuploidy.score)

# get 8p neutral cell lines
df_aneu_score %>% 
  filter(X8p == 0) %>% 
  dplyr::select(X) %>% 
  unlist(use.names = F) -> chr_8p_disomy

# median aneuploidy for chr8p neutral
aneuploid_dep %>% 
  dplyr::filter(DepMap_ID  %in% chr_8p_disomy) -> chr8p_neutral_flt
median(chr8p_neutral_flt$Aneuploidy.score)

# checking for the aneuploidy score for the cell lines
# aneuploidy score distribution
aneuploid_dep %>% 
  filter(DepMap_ID %in% c(chr_8p_loss,chr_8p_disomy )) %>% 
  mutate(aneu_con = ifelse(DepMap_ID %in% chr_8p_loss, "monosomy", "disomy")) %>% 
  ggplot(aes(x = Aneuploidy.score)) +
  geom_vline(xintercept = 7,linetype = "dashed") +
  geom_vline(xintercept = 27,linetype = "dashed") +
  geom_rangeframe()+
  theme_tufte()+
  geom_density(aes(fill = aneu_con,alpha = 0.5))+
  scale_fill_manual(values = c( colors[2],colors[1]))

# whole genome duplication
# the diffence of whole genome doubling in 1q trisomy and 1q disomy
aneuploid_dep %>% 
  filter(DepMap_ID %in% c(chr_8p_loss,chr_8p_disomy )) %>% 
  mutate(aneu_con = ifelse(DepMap_ID %in% chr_8p_loss, "monosomy", "disomy")) %>% 
  group_by(aneu_con) %>% 
  mutate(sum_group = n()) %>% ungroup() %>% 
  group_by(Genome.doublings,aneu_con) %>% 
  mutate(number = n()) %>% 
  distinct(Genome.doublings, .keep_all = T) %>% 
  mutate(ratio = number/ sum_group) %>% 
  ungroup() %>% 
  ggplot(aes(x = Genome.doublings, y=ratio, fill = aneu_con)) +
  geom_histogram(position = "dodge",stat = "identity")+
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c( colors[2],colors[1]))


# before the cut off and check how the deviation between 1q trisomy and 1q disomy
df_aneu_score %>% 
  filter(X %in% c(chr_8p_loss,chr_8p_disomy )) %>% 
  mutate(aneu_con = ifelse(X %in% chr_8p_loss, "monosomy", "disomy"))-> df_aneu_score_bef_flt
correct_aneu_bef_flt = lapply(c(1:(dim(df_aneu_score_bef_flt)[2]-2)),function(x){
  df_col = df_aneu_score_bef_flt[,x+1]
  df_col_corr = df_col[df_col == -1] = 1
  col_name = colnames(df_aneu_score_bef_flt)[x +1]
  df = data.frame(df_col)
  colnames(df) = col_name
  return(df)
})

# combine the dataset
combin_df_bef = do.call(cbind,correct_aneu_bef_flt)
combin_df_bef$DepID = df_aneu_score_bef_flt$X
combin_df_bef$aneu_con =  df_aneu_score_bef_flt$aneu_con

# plot the figures 
combin_df_bef %>% 
  group_by(aneu_con) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  gather(X1p:X22q,key = "chr_arm", value = "value") %>% 
  #filter(chr_arm != "X1q") %>% 
  
  group_by(aneu_con, chr_arm) %>% 
  
  mutate(aneuploidy_sum = sum(value)) %>% 
  distinct(chr_arm, .keep_all = T) %>% 
  mutate(ratio = aneuploidy_sum/ total) %>% 
  ungroup() %>% 
  mutate(chr_arm = factor(chr_arm, levels = colnames(combin_df_bef)[1:39])) %>% 
  ggplot(aes(x = chr_arm, y = ratio, shape =aneu_con, color =aneu_con, group=aneu_con )) +
  #geom_histogram(position = "dodge", stat = "identity") +
  geom_point() + geom_smooth(span = .2)+
  
  theme_classic()+
  scale_color_manual(values = c( colors[2],colors[1]))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))-> p_bef




# based on the density plot, we found that :
# 1. to make the accurate comparison, I need first to guarantee my cells similar genomic background 
# 2. If I compare the cell with different genomic background, this could drive a false or unrelevant target
# 3. because of our isogeneic cells, with high purity and less interferece, the better
# 4.  set up the cut off at aneuploidy > 7 - 27 as our primary setting
# 5. filter out whole genome doubling cells

# Primary setting
# based on the aneuploidy score
# 8p monosomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_8p_loss) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 27) %>%
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_8p_loss_pri_flt

# 1q disomy: after primary filter
aneuploid_dep %>% 
  filter(DepMap_ID %in% chr_8p_disomy) %>% 
  filter(Aneuploidy.score > 7 & Aneuploidy.score < 27) %>% 
  filter(Genome.doublings != 2) %>% 
  dplyr::select(DepMap_ID) %>% 
  unlist(use.names = F) -> chr_8p_diso_pri_flt  

# plot after filtreation
aneuploid_dep %>% 
  filter(DepMap_ID %in% c(chr_8p_loss,chr_8p_disomy )) %>% 
  mutate(aneu_con = ifelse(DepMap_ID %in% chr_8p_loss, "monosomy", "disomy")) %>% 
  dplyr::filter(DepMap_ID %in% c(chr_8p_loss_pri_flt, chr_8p_diso_pri_flt)) %>% 
  ggplot(aes(x = Aneuploidy.score)) +
  #geom_vline(xintercept = 7,linetype = "dashed") +
  #geom_vline(xintercept = 27,linetype = "dashed") +
  geom_rangeframe()+
  theme_tufte()+
  geom_density(aes(fill = aneu_con,alpha = 0.5))+
  scale_fill_manual(values = c( colors[2],colors[1]))

# after primary filter, I need to check how geneal aneuploidy outside chromosome 18q looks like
# I first need to define the aneuploidy condition in different chromosome 
# -1 or 1 as aneuploidy all equal to 1
# now we have to check if the cell have similarity on the other chromosome gain or loss
df_aneu_score %>% 
  filter(X %in% c(chr_8p_loss_pri_flt,chr_8p_diso_pri_flt )) %>% 
  mutate(aneu_con = ifelse(X %in% chr_8p_loss_pri_flt, "monosomy", "disomy"))-> df_aneu_score_pri_flt
correct_aneu_pri_flt = lapply(c(1:(dim(df_aneu_score_pri_flt)[2]-2)),function(x){
  df_col = df_aneu_score_pri_flt[,x+1]
  df_col_corr = df_col[df_col == -1] = 1
  col_name = colnames(df_aneu_score_pri_flt)[x +1]
  df = data.frame(df_col)
  colnames(df) = col_name
  return(df)
})

# combine the dataset
combin_df = do.call(cbind,correct_aneu_pri_flt)
combin_df$DepID = df_aneu_score_pri_flt$X
combin_df$aneu_con =  df_aneu_score_pri_flt$aneu_con


# need to have check if there is difference of fraction between 7p gain and 7p loss  
#View(df_aneu_score_pri_flt)
fraction_gain_loss = lapply(c(1:39), function(x){
  df_1 = df_aneu_score_pri_flt[,x + 1]
  # trisomy
  tri = df_1[df_aneu_score_pri_flt$aneu_con == "monosomy"]
  tri_gain = sum(tri == 1) / length(tri)
  tri_loss = sum(tri == -1)/ length(tri)
  
  # disomy
  diso = df_1[df_aneu_score_pri_flt$aneu_con == "disomy"]
  di_gain = sum(diso == 1) / length(diso)
  di_loss = sum(diso == -1)/ length(diso)
  
  frac = c(tri_gain,tri_loss, di_gain, di_loss)
  aneu_con = c(rep(c("gain","loss"),time = 2))
  chr = colnames(df_aneu_score_pri_flt)[x + 1]
  tri_di = c(rep(c("monosomy","disomy"),each = 2))
  
  df_2 = data.frame(frac,aneu_con,chr,tri_di)
  
  return(df_2)
  
  
})

fraction_gain_loss_comb = do.call(rbind,fraction_gain_loss)
colnames(fraction_gain_loss_comb) = c("value", "aneu_con", "chr","tri_di")

# plot the fraction changes
# I am not sure if i need extra weight system 
# looks i do not need it 
chr_name = colnames(df_aneu_score_pri_flt)[2:40]
order_chr = c()
for(i in c(1:39)){
  order_chr[i] = chr_name[39-i + 1]
}

# plot the distribution difference based on the chromosome gain or loss
fraction_gain_loss_comb %>% 
  mutate(chr = factor(chr, levels = order_chr))->p
p %>% 
  ggplot(aes(x = chr) ) +
  geom_bar(data = subset( p,tri_di =="disomy"), aes(y = value, fill = aneu_con), 
           stat = "identity",position = "stack", alpha = .9)  +
  geom_bar(data = subset( p,tri_di =="monosomy"), aes(y = -value, fill = aneu_con), 
           stat = "identity",position = "stack", alpha = .9)+
  #facet_wrap(~tri_di) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  coord_flip() +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = colors)
#facet_grid(~tri_di)


# plot the figures 
combin_df %>% 
  group_by(aneu_con) %>% 
  mutate(total = n()) %>% 
  ungroup() %>% 
  gather(X1p:X22q,key = "chr_arm", value = "value") %>% 
  #filter(chr_arm != "X1q") %>% 
  
  group_by(aneu_con, chr_arm) %>% 
  
  mutate(aneuploidy_sum = sum(value)) %>% 
  distinct(chr_arm, .keep_all = T) %>% 
  mutate(ratio = aneuploidy_sum/ total) %>% 
  ungroup() %>% 
  mutate(chr_arm = factor(chr_arm, levels = colnames(combin_df_bef)[1:39])) %>% 
  ggplot(aes(x = chr_arm, y = ratio, shape =aneu_con, color =aneu_con, group=aneu_con )) +
  #geom_histogram(position = "dodge", stat = "identity") +
  geom_point() + geom_smooth(span = .2)+
  theme_classic()+
  scale_color_manual(values = c( colors[2],colors[1]))+
  #scale_fill_manual(values = c( colors[2],colors[1]))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) ->p_aft

# combine the figure
ggarrange(p_bef,p_aft ,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)




# 8p loss dependendcy
# Gene dependency 24Q2 data set
df_gene_dep = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /CRISPRGeneEffect.csv")

# cell line info
metasheet_cellline = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /Model.csv")

# change the name of data
colnames(df_gene_dep)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_Dep = sub("\\..*","",colnames(df_gene_dep)[2:dim(df_gene_dep)[2]])

# replace gene to colnames
colnames(df_gene_dep)[2:dim(df_gene_dep)[2]] = genelist_Dep
head(df_gene_dep[,1:6])

#(dim(df_gene_dep)[2]-1)

chr_18q_NormLRT = lapply(c(1:(dim(df_gene_dep)[2]-1)), function(x){
  
  depscore = df_gene_dep[, c(1, x + 1)]
  gene = colnames(depscore)[2]
  chr18_q_depscore = depscore[depscore$DepMap_ID %in% chr_8p_loss_pri_flt,]
  chr18_q_nondepscore = depscore[depscore$DepMap_ID %in% chr_8p_diso_pri_flt,]
  chr18q_median = median(chr18_q_depscore[,2][!is.na(chr18_q_depscore[,2])])
  chr18q_non_median = median(chr18_q_nondepscore[,2][!is.na(chr18_q_nondepscore[,2])])
  diff_depscore =  chr18q_non_median - chr18q_median
  #print(x)
  # calculate NormLRT
  ##load function
  #Log-likelihood for the skew normal
  llsn <- function(par){
    if(par[2]>0) return( sum(dsn(data.sim, par[1], par[2], par[3], log = T))  )
    else return(-Inf)
  }
  
  data.sim = depscore[,2][!is.na(depscore[,2])]
  init1 <- c(mean(data.sim), sd(data.sim),-0.5)
  init2 <- c(mean(data.sim), sd(data.sim),0)
  init3 <- c(mean(data.sim), sd(data.sim),0.5)
  
  OPT1 <- optim(init1, llsn, control = list(fnscale = -1, maxit = 10000))
  OPT2 <- optim(init2, llsn, control = list(fnscale = -1, maxit = 10000))
  OPT3 <- optim(init3, llsn, control = list(fnscale = -1, maxit = 10000))
  
  indmax <- which.max(c(OPT1$value,OPT2$value,OPT3$value))
  indmax # Best ini
  if(indmax==1) MLESN <- OPT1$par; if(indmax==2) MLESN <- OPT2$par; if(indmax==3) MLESN <- OPT3$par;
  
  
  #Log-likelihood for the normal 
  lln <- function(par){
    if(par[2]>0) return( sum(dnorm(data.sim, par[1], par[2], log = T))  )
    else return(-Inf)
  }
  
  # MLE for the normal parameters
  MLEN <- optim(c(0,1), lln, control = list(fnscale = -1))$par
  
  LRT <- -2*(lln(MLEN) - llsn(MLESN))
  
  #make a table
  data = data.frame(gene, diff_depscore, LRT, chr1_median, chr1_non_median)
})

# integrate the list
NormLRT_chr7p = do.call(rbind, chr_7p_NormLRT)
colnames(NormLRT_chr1q) = c("gene", "diff", "NormLRT", "chr1_median","chr1_non_median")
#write_xlsx(NormLRT_diff_1,path = "~/Documents/neuroblastoma database/DEP_MAP data/DEPmap cripsr analysis final/depNormLRT.xlsx")
write_xlsx(NormLRT_chr7p, "chr7p_dep_screen.xlsx")


# compare the diff wilcox
chr_8p_wilcox = lapply(c(1:(dim(df_gene_dep)[2]-1)), function(x){
  
  depscore = df_gene_dep[, c(1, x + 1)]
  gene = colnames(depscore)[2]
  chr8p_depscore = depscore[depscore$DepMap_ID %in% chr_8p_loss_pri_flt,]
  chr8p_nondepscore = depscore[depscore$DepMap_ID %in% chr_8p_diso_pri_flt,]
  chr8p_median = median(chr8p_depscore[,2][!is.na(chr8p_depscore[,2])])
  chr8p_non_median = median(chr8p_nondepscore[,2][!is.na(chr8p_nondepscore[,2])])
  diff_depscore =  chr8p_non_median - chr8p_median
  #print(x)
  # calculate NormLRT
  if(sum(!is.na(chr8p_depscore[,2])) >2 & sum(!is.na(chr8p_nondepscore[,2])) >2){
    p_value = wilcox.test(chr8p_depscore[,2], chr8p_nondepscore[,2],exact = FALSE, correct = FALSE)$p.value
    
  }else{
    p_value = 1
    
  }
  
  
  #make a table
  data = data.frame(gene, diff_depscore, p_value, chr8p_median, chr8p_non_median)
})

#combine the data
wilcox_chr8p_loss = do.call(rbind, chr_8p_wilcox)
colnames(wilcox_chr8p_loss) = c("gene", "diff", "p_value", "chr8p_median","chr8p_non_median")

#save the data
wilcox_chr8p_loss$chr = "chr8_p"

#path = "~/Documents/Postdoc_sheltzer/Publication/Paralog dependency/dependency_summary/all datasheet/"

#write.csv(wilcox_chr8p_loss, file.path(path,"chr8_p_dep.csv" ), row.names = F)


wilcox_chr8p_loss %>% 
  dplyr::filter(diff >0.1) %>% 
  dplyr::filter(p_value < 0.01 & chr8p_median < -0.3) -> sig_selected_gene

# PPP2CA and PPP2CB (chr8p)

sig_selected_gene$chr = "chr8_p"

path = "~/Documents/Postdoc_sheltzer/Publication/Paralog dependency/dependency_summary/"

write.csv(sig_selected_gene, file.path(path,"chr8_p_dep.csv" ), row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)

table1<- bitr(sig_selected_gene$gene, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = table1$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)

View(ego[1:200,])



# Mahattan plot
library(biomaRt)
library(Homo.sapiens)

#retrieve the gene start postion
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

results_gene_dep <- getBM(attributes = c("chromosome_name",  "hgnc_symbol","start_position"),
                          filters = "hgnc_symbol", values = genelist_Dep, mart = mart)

# get the gene on the chromosme 1:22 and X,Y
chr_name = c(1:22,"X","Y")
results_gene_dep %>% 
  dplyr::filter(chromosome_name %in% chr_name) %>% 
  dplyr::rename(gene = hgnc_symbol)-> gene_chrname
# there are some genes does not include into the study
# get the chromosome length

#chr_length = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)

#chr_length_flt = unname(chr_length[1:24])

# get the chr length
cytoband_chr = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Dataset/cytoBand.txt",
                          header = F)

# rename the file
colnames(cytoband_chr) = c("chr","start","end","band", "other")

# get the chromosome info
chr_name_2 = paste("chr", c(1:22,"X","Y"), sep = "")

# get the chromsome length table
start_list = c()
end_list = c()
length_list = c()


for (i in chr_name_2) {
  df = cytoband_chr[cytoband_chr$chr == i,]
  # start
  start = min(df$start)
  start_list = c(start_list, start)
  # end
  end = max(df$end)
  end_list = c(end_list, end)
  # length chr
  length_chr = end - start
  length_list = c(length_list, length_chr)
  
  
}
chr_length = length_list
chr_length_flt = length_list

# calculate acculated chromosome
acu_length = c()
for(i in c(1:24)){
  if(length(chr_length_flt[i-1]) == 0 ){
    length_chr = 0
    acu_length = c(acu_length , length_chr)
  }else{
    length_chr = chr_length_flt[i-1] + acu_length[i-1]
    acu_length = c(acu_length , length_chr)
  }
}

chr_query = data.frame(chr_name, acu_length)
colnames(chr_query) = c("chromosome_name", "chr_length")


# color annotation
colors = rep(c("blue","grey"), time = 12)
chr_color = data.frame(chr_name, colors)
colnames(chr_color) = c("chromosome_name", "color1")

# combine the data
levels_chr = c(1:22,"X","Y")
wilcox_chr8p_loss %>% 
  inner_join(gene_chrname, by = "gene") %>% 
  left_join(chr_query,  by = "chromosome_name") %>% 
  left_join(chr_color, by = "chromosome_name") %>% 
  mutate(bp_cum = chr_length + start_position) -> p1

#library
library(ggthemes)
library(ggpubr)

# bp_cum cell
p1 %>%   
  group_by(chromosome_name) %>% 
  summarize(center = mean(bp_cum)) %>% 
  ungroup() -> p2

# set the chromosme center
p1 %>% 
  left_join(p2, by = "chromosome_name") ->p


# library
library(ggrepel)
# color set
color_set= c("#EF8536","#3A76AF")

p %>% 
  mutate(color2 = ifelse(p_value < 0.01 & diff > 0.1 &  chr8p_median < -0.3 ,"orange", color1)) %>% 
  mutate(color = ifelse(p_value < 0.01 & diff < 0, "black", color2)) %>% 
  mutate(color = factor(color,levels = c("orange","black", "blue","grey"))) %>% 
  arrange(desc(color)) %>% 
  dplyr::filter(diff > 0) %>% 
  ggplot(aes(x = bp_cum, y = -log10(p_value),
             size = -log10(p_value),
             label = gene))+
  geom_point(aes(color = color),
             alpha = 0.75) +
  geom_text_repel(
    data = subset(p, gene %in% c("PPP2CA","POLR3D","ATP6V1H")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
    
  )+
  scale_x_continuous(
    label = p$chromosome_name,
    breaks = p$center
  ) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,9)) +
  scale_color_manual(values = c("#EF8536" ,"#D7BFDC","#7852A9")) +
  scale_size_continuous(range = c(0.5, 3)) +
  geom_rangeframe()+
  theme_tufte()+
  theme(legend.position  = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x = "chromosome", 
       y = "-log10(p value)")


# correlation dependency
## analysis here: correlation and dep and exp
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

# PPP2CA is dependent on chr18q loss

PPP2CA_dep = df_gene_dep[,c(1,which(colnames(df_gene_dep) == "PPP2CA"))]
PPP2CA_dep %>% 
  dplyr::filter(DepMap_ID %in% c(chr_8p_diso_pri_flt,chr_8p_loss_pri_flt)) -> PPP2CA_dep_flt

PPP2CA_dep_flt %>% 
  dplyr::mutate(condition = ifelse(DepMap_ID %in%chr_8p_loss_pri_flt ,"8p_loss", "8p_non_loss")) -> PPP2CA_dep

# plot the dependency 
PPP2CA_dep %>% 
  ggplot(aes(x = condition, y = PPP2CA)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(color_set[1],"#7852A9"))+
  scale_color_manual(values = c(color_set[1],"#7852A9")) +
  theme(legend.position = "none") +
  labs(y = "PPP2CA_dep") -> P

P
# p value = 3.414e-09
wilcox.test(PPP2CA_dep[PPP2CA_dep$condition == "8p_loss",]$PPP2CA, PPP2CA_dep[PPP2CA_dep$condition == "8p_non_loss",]$PPP2CA)


# dependency and exp correlation analysis
PPP2CA_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  PPP2CA_dep_flt %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,2],comb_tmp[,3])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})

# combine the list
PPP2CA_exp_corr_cmb = do.call( rbind,PPP2CA_exp_corr)
colnames(PPP2CA_exp_corr_cmb) = c("gene", "p_value", "corr")


library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
PPP2CA_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(PPP2CA_exp_corr_cmb,gene%in% c("PPP2CB")),
    aes(label = gene),
    size = 4,
    box.padding = 0.6,
    point.padding = 0.5,
    #label.padding = 0.25,
    min.segment.length = unit(0, 'lines'),
    max.overlaps = 20,
    col="black"
  ) +
  geom_rangeframe()+
  theme_tufte() +
  scale_color_manual(values = c(colors[2],"grey",colors[1])) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 15),
         legend.position = "none")



# expression and dependency correlation

PPP2CB_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "PPP2CB"))]
PPP2CA_dep %>% 
  left_join(PPP2CB_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(PPP2CA)) %>% 
  filter(!is.na(PPP2CB)) -> cmb_df
cmb_df$pc = predict(prcomp(~PPP2CA+PPP2CB, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = PPP2CB, y = PPP2CA)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  #scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  scale_color_gradient2(low = "#FCAE1E" ,mid = "#E2DDFC",  high = "#543D7B") + 
  theme(legend.position = "none") +
  labs( y = "PPP2CA_dep (CRISPR)", x = "PPP2CB exp (TPM)") +
  geom_text( x = 10, y = -0.5, label = "R = 0.2864524")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 6.799e-06", x = 10, y =-0.55) 

# correlation tesPPP2CB# correlation test
#p-value = 2.01e-08
cor_result = cor.test(cmb_df$PPP2CA, cmb_df$PPP2CB)

# PPP2CB expression difference 

PPP2CB_exp %>% 
  dplyr::filter(DepMap_ID %in% c(chr_8p_diso_pri_flt,chr_8p_loss_pri_flt)) -> PPP2CB_exp_flt

PPP2CB_exp_flt %>% 
  dplyr::mutate(condition = ifelse(DepMap_ID %in%chr_8p_loss_pri_flt ,"8p_loss", "8p_non_loss")) -> PPP2CB_exp_anno

# plot the dependency 
PPP2CB_exp_anno %>% 
  ggplot(aes(x = condition, y = PPP2CB)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(color_set[1],"#7852A9"))+
  scale_color_manual(values = c(color_set[1],"#7852A9")) +
  theme(legend.position = "none") +
  labs(y = "PPP2CB expression") -> P2

P2


# p value = 3.414e-09
t.test(PPP2CB_exp_anno[PPP2CB_exp_anno$condition == "8p_loss",]$PPP2CB, PPP2CB_exp_anno[PPP2CB_exp_anno$condition == "8p_non_loss",]$PPP2CB)

# expression difference analysis
expression_diff_ana = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  #print(x)
  df = CCLE_exp[,c(1, x + 1)]
  df %>% 
    dplyr::filter(DepMap_ID %in% c(chr_18q_diso_pri_flt,chr_18q_loss_sec_flt)) %>% 
    dplyr::mutate(aneu_com = ifelse(DepMap_ID %in% chr_18q_diso_pri_flt, 
                                    "18q_disomy", "18q_monosomy" )) ->P
  gene = colnames(CCLE_exp)[x +1]
  #get the median
  median_chr18q_loss = median(P[P$aneu_com == "18q_monosomy",][,2][!is.na(P[P$aneu_com == "18q_monosomy",][,2])])
  median_chr18q_diso = median(P[P$aneu_com == "18q_disomy",][,2][!is.na(P[P$aneu_com == "18q_disomy",][,2])])
  if(sum(!is.na(P[P$aneu_com == "18q_disomy",][,2])) > 1 & sum(!is.na(P[P$aneu_com == "18q_monosomy",][,2])) > 1 ){
    p_value = t.test(P[P$aneu_com == "18q_disomy",][,2], P[P$aneu_com == "18q_monosomy",][,2])$p.value
    
  }else{
    p_value = 1
  }
  com_df = data.frame(gene, p_value,median_chr18q_loss,median_chr18q_diso)
  return(com_df)
  
})
com_expression_diff_ana = do.call(rbind, expression_diff_ana)
colnames(com_expression_diff_ana) = c("gene", "p_value", "median_chr18q_loss", "median_chr18q_diso")

# check the psotive corr gene
com_expression_diff_ana %>% 
  dplyr::filter(median_chr18q_loss < median_chr18q_diso) %>% 
  dplyr::filter(p_value < 0.01) %>%
  dplyr::select(gene) %>% unlist(use.names = F) ->pos_sig_gene

PIK3C3_exp_corr_cmb %>% 
  dplyr::filter(corr >0.15,p_value < 0.01) %>% 
  dplyr::select(gene) %>% unlist(use.names = F) ->pos_sig_corr

pos_cor_genelist = pos_sig_gene[pos_sig_gene%in%pos_sig_corr]

library(clusterProfiler)
library(org.Hs.eg.db)

table1<- bitr(pos_cor_genelist, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = table1$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)

View(ego[1:200,])


# check the negativ corr gene
com_expression_diff_ana %>% 
  dplyr::filter(median_chr18q_loss > median_chr18q_diso) %>% 
  dplyr::filter(p_value < 0.01) %>%
  dplyr::select(gene) %>% unlist(use.names = F) ->neg_sig_gene

PIK3C3_exp_corr_cmb %>% 
  dplyr::filter(corr < -0.15,p_value < 0.01) %>% 
  dplyr::select(gene) %>% unlist(use.names = F) ->neg_sig_corr

neg_cor_genelist = neg_sig_gene[neg_sig_gene%in%neg_sig_corr]

# function analysis
table1<- bitr(neg_cor_genelist, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = table1$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

View(ego[1:200,])


#EXPRESSION DIFFERENCE of DDX3Y in chrY
PPP2CB_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) %in% c("PPP2CB")))]
PPP2CB_exp %>% 
  dplyr::filter(DepMap_ID %in% c(chr_8p_diso_pri_flt,chr_8p_loss_pri_flt)) %>% 
  dplyr::mutate(aneu_com = ifelse(DepMap_ID %in% chr_8p_loss_pri_flt, "chr8p_monosomy", "chr8p_disomy" )) ->P
#left_join(cell_info_flt, by ="DepMap_ID" ) %>% 
#dplyr::filter(Gender %in% c("male", "female")) ->p
P %>% 
  ggplot(aes(x = factor(aneu_com, levels = c("chr8p_monosomy", "chr8p_disomy")), y = PPP2CB )) +
  geom_boxplot(aes(color = aneu_com,fill = aneu_com,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = aneu_com),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c("#7852A9", colors[1]))+
  scale_color_manual(values = c("#7852A9", colors[1])) +
  theme(legend.position = "none") +
  labs(x = "", y = "PPP2CB expression") -> p2

#ggarrange( p2, ncol = 2,nrow = 1)
p2
# p = 3.263124e-19
p_value = t.test(P[P$aneu_com == "chr8p_monosomy",]$PPP2CB, P[P$aneu_com == "chr8p_disomy",]$PPP2CB)
p_value$p.value


# copy number of PPP2CB
# load the cnv data
CNV_gene = read.csv("~/Documents/Postdoc_sheltzer/Depmp data /PortalOmicsCNGeneLog2.csv")

# change the name of CCLE CNV data
colnames(CNV_gene)[1] = "DepMap_ID"

# gene the gene name 
genelist_cnv = sub("\\..*","",colnames(CNV_gene)[2:dim(CNV_gene)[2]])

# replace gene to colnames
colnames(CNV_gene)[2:dim(CNV_gene)[2]] = genelist_cnv

# get PPP2CB cnv
PPP2CB_cnv = CNV_gene[,c(1, which(colnames(CNV_gene) %in% c("PPP2CB")))]
colnames(PPP2CB_cnv)[2] = "PPP2CB_cnv"
P %>% 
  left_join(PPP2CB_cnv, by = "DepMap_ID") ->PPP2CB_merge
PPP2CB_merge %>% 
  ggplot(aes(x = factor(aneu_com, levels = c("chr8p_monosomy", "chr8p_disomy")), y = PPP2CB_cnv )) +
  geom_boxplot(aes(color = aneu_com,fill = aneu_com,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = aneu_com),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  ylim(0,1.5)+
  theme(legend.position = "none") +
  labs(x = "", y = "PPP2CB copy number") -> p2

# arrange the plot
ggarrange(p2, ncol = 2,nrow = 1)


# copy number difference
# p_value = 1.430381e-53
p_value = t.test(PPP2CB_merge[PPP2CB_merge$aneu_com == "chr8p_monosomy",]$PPP2CB_cnv, PPP2CB_merge[PPP2CB_merge$aneu_com == "chr8p_disomy",]$PPP2CB_cnv)
p_value$p.value

# correlation between copy number and expression
PPP2CB_merge %>% 
  dplyr::filter(!is.na(PPP2CB_cnv)) ->PPP2CB_merge_flt

#weigth
PPP2CB_merge_flt$pc = predict(prcomp(~PPP2CB_cnv+PPP2CB, PPP2CB_merge_flt))[,1]

#cor test
cor.test(PPP2CB_merge_flt$PPP2CB_cnv,PPP2CB_merge_flt$PPP2CB)


PPP2CB_merge_flt %>% 
  ggplot(aes(x = PPP2CB_cnv, y = PPP2CB)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2)+
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  labs(y = "PPP2CB expression", x = "PPP2CB copy number")+
  scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = 1)) +
  theme(legend.position = "none") +
  geom_text( x = 1.8, y =4 , label = "R = 0.5730618 ")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value  < 2.2e-16", x = 1.8, y =3.8)




