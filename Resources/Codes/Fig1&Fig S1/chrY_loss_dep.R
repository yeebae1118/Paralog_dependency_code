# project : chr Y loss dependency
# Date: Oct.30.2024
# Author: Yi

# library
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(org.Hs.eg.db)
library(biomaRt)


# data loading
# gene copy number
CCLE_gene_cnv = read.csv("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /PortalOmicsCNGeneLog2.csv")

# cell and patient info
cell_line_id = read.delim("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Depmp data /aneuploid depmap/Cell_lines_annotations_20181226.txt")

# rename the gene name on CCLE_gene_CNV
# change the name of data
colnames(CCLE_gene_cnv)[1] = "DepMap_ID"
#dim(df_gene_dep)

# gene the gene name 
dim(CCLE_gene_cnv)
genelist_Dep = sub("\\..*","",colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]])

# replace gene to colnames
colnames(CCLE_gene_cnv)[2:dim(CCLE_gene_cnv)[2]] = genelist_Dep
head(CCLE_gene_cnv[,1:6])

# retrieve Y chromosome genes
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
results <- getBM(attributes = c("chromosome_name", "entrezgene_id", "hgnc_symbol"),
                 filters = "chromosome_name", values = "Y", mart = mart)

Y_chr_genename = results$hgnc_symbol[results$hgnc_symbol != ""]

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


# Y loss dependendcy
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

# compare the diff wilcox
chr_chrY_wilcox = lapply(c(1:(dim(df_gene_dep)[2]-1)), function(x){
  
  depscore = df_gene_dep[, c(1, x + 1)]
  gene = colnames(depscore)[2]
  chrY_depscore = depscore[depscore$DepMap_ID %in% Y_loss_sample,]
  chrY_nondepscore = depscore[depscore$DepMap_ID %in% Y_normal,]
  chrY_median = median(chrY_depscore[,2][!is.na(chrY_depscore[,2])])
  chrY_non_median = median(chrY_nondepscore[,2][!is.na(chrY_nondepscore[,2])])
  diff_depscore =  chrY_non_median - chrY_median
  #print(x)
  # calculate NormLRT
  if(sum(!is.na(chrY_depscore[,2])) >2 & sum(!is.na(chrY_nondepscore[,2])) >2){
    p_value = wilcox.test(chrY_depscore[,2], chrY_nondepscore[,2],exact = FALSE, correct = FALSE)$p.value
    
  }else{
    p_value = 1
    
  }
  
  
  #make a table
  data = data.frame(gene, diff_depscore, p_value, chrY_median, chrY_non_median)
})

#combine the data
wilcox_chrY_loss = do.call(rbind, chr_chrY_wilcox)
colnames(wilcox_chrY_loss) = c("gene", "diff", "p_value", "chrY_median","chrY_non_median")

# save the data
wilcox_chrY_loss$chr = "chrY"

path = "~/Documents/Postdoc_sheltzer/Publication/Paralog dependency/dependency_summary/all datasheet/"

write.csv(wilcox_chrY_loss, file.path(path,"chrY_dep.csv" ),row.names = F )


wilcox_chrY_loss %>% 
  dplyr::filter(diff >0.1) %>% 
  dplyr::filter(p_value < 0.01 & chrY_median < -0.3) ->  sig_selected_gene
sig_selected_gene$chr = "chrY"

path = "~/Documents/Postdoc_sheltzer/Publication/Paralog dependency/dependency_summary/"

write.csv(sig_selected_gene, file.path(path,"chrY_dep.csv" ),row.names = F )

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
wilcox_chrY_loss %>% 
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
  mutate(color2 = ifelse(p_value < 0.01 & diff > 0.15 &chrY_median < -0.3 ,"orange", color1)) %>% 
  mutate(color = ifelse(p_value < 0.01 & diff < -0.15, "black", color2)) %>% 
  mutate(color = factor(color,levels = c("orange","black", "blue","grey"))) %>% 
  arrange(desc(color)) %>% 
  dplyr::filter(diff > 0) %>% 
  ggplot(aes(x = bp_cum, y = -log10(p_value),
             size = -log10(p_value),
             label = gene))+
  geom_point(aes(color = color),alpha = 0.75) +
  geom_text_repel(
    data = subset(p, gene %in% c("DDX3X","EIF1AX","RPS4X")),
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 16)) +
  scale_color_manual(values = c("#EF8536" ,"#D7BFDC","#7852A9")) +
  scale_size_continuous(range = c(0.5, 3)) +
  geom_rangeframe()+
  theme_tufte()+
  theme(legend.position  = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x = "chromosome", 
       y =  "-log10(p value)")


# DDX3X is dependent on Y loss

DDX3X_dep = df_gene_dep[,c(1,which(colnames(df_gene_dep) == "DDX3X"))]
DDX3X_dep %>% 
  dplyr::filter(DepMap_ID %in% c(Y_loss_sample,Y_normal)) -> DDX3X_dep_flt

DDX3X_dep_flt %>% 
  dplyr::mutate(condition = ifelse(DepMap_ID %in%Y_loss_sample ,"Y_loss", "Y_non_loss")) -> DDX3X_dep

# plot the dependency 
DDX3X_dep %>% 
  ggplot(aes(x = condition, y = DDX3X)) +
  geom_boxplot(aes(color = condition,fill = condition,alpha = 0.3),width = 0.25,outlier.shape = NA) +
  geom_jitter( aes(color = condition),position=position_jitter(0.12),size = 1,alpha = 0.3) +
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],"#7852A9"))+
  scale_color_manual(values = c(colors[1],"#7852A9")) +
  theme(legend.position = "none") +
  labs(y = "DDX3X_dep") -> P

P

# p value = 3.414e-09
wilcox.test(DDX3X_dep[DDX3X_dep$condition == "Y_loss",]$DDX3X, DDX3X_dep[DDX3X_dep$condition == "Y_non_loss",]$DDX3X)

# correlation dependency
## analysis here: correlation and dep and exp
# retrieve gene expression data
CCLE_exp = read.csv("~/Documents/Postdoc_sheltzer/Depmp data /OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")

# change the name of CCL2_exp data
colnames(CCLE_exp)[1] = "DepMap_ID"
dim(df_gene_dep)

# gene the gene name 
dim(df_gene_dep)
genelist_exp = sub("\\..*","",colnames(CCLE_exp)[2:dim(CCLE_exp)[2]])

# replace gene to colnames
colnames(CCLE_exp)[2:dim(CCLE_exp)[2]] = genelist_exp

# retrieve DDX3X dependency in 7p gain cancer

DDX3X_dep = df_gene_dep[,c(1,which(colnames(df_gene_dep) == "DDX3X"))]
DDX3X_dep %>% 
  dplyr::filter(DepMap_ID %in% c(Y_loss_sample,Y_normal)) -> DDX3X_dep_flt

# dependency and exp correlation analysis
DDX3X_exp_corr = lapply(c(1:(dim(CCLE_exp)[2]-1)), function(x){
  exp_df = CCLE_exp[,c(1,x +1)]
  DDX3X_dep_flt %>% 
    left_join(exp_df, by = "DepMap_ID") -> comb_tmp
  statistic_ana = cor.test(comb_tmp[,2],comb_tmp[,3])
  p_value = statistic_ana$p.value
  cor = statistic_ana$estimate
  gene = colnames(exp_df)[2]
  df_retn = data.frame(gene, p_value,cor)
  return(df_retn)
})

# combine the list
DDX3X_exp_corr_cmb = do.call( rbind,DDX3X_exp_corr)
colnames(DDX3X_exp_corr_cmb) = c("gene", "p_value", "corr")




library(ggrepel)
# color set
colors = c("#EF8536","#3A76AF")
DDX3X_exp_corr_cmb %>% 
  dplyr::mutate(color_order_1 = ifelse(p_value > 0.05,"grey","orange")) %>% 
  dplyr::mutate(color_order_2 = ifelse(p_value < 0.05 & corr < 0,"blue",color_order_1)) %>% 
  dplyr::mutate(color_order_3 = ifelse(p_value < 0.05 & corr >0 ,"orange",color_order_2)) %>% 
  ggplot(aes(x = corr, y= -log10(p_value), size =  -log10(p_value),alpha = abs(corr)))+
  geom_point(aes(color = color_order_3))+
  geom_text_repel(
    data = subset(DDX3X_exp_corr_cmb,gene%in% c("DDX3Y","UTY", "EIF1AY", "KDM5D","PRKY",
                                                "USP9Y","RPS4Y1","ZFY","NLGN4Y","TMSB4Y")),
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

DDX3Y_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) == "DDX3Y"))]
DDX3X_dep%>% 
  left_join(DDX3Y_exp, by = "DepMap_ID") -> cmb_df


#weigth
cmb_df %>% 
  filter(!is.na(DDX3X)) %>% 
  filter(!is.na(DDX3Y)) -> cmb_df
cmb_df$pc = predict(prcomp(~DDX3X+DDX3Y, cmb_df))[,1]

cmb_df %>% 
  ggplot(aes(x = DDX3Y, y = DDX3X)) +
  geom_point(aes(color = pc),alpha = 0.5,size = 2) +
  geom_smooth(method = "lm") +
  geom_rangeframe()+
  theme_tufte() +
  #scale_color_gradientn(colours = met.brewer(name = "Hokusai3", n= 7,type = "continuous",direction = -1)) +
  scale_color_gradient2(high = "#FCAE1E" ,mid = "#E2DDFC",  low = "#543D7B") + 
  theme(legend.position = "none") +
  labs( y = "DDX3X_dep (CRISPR)", x = "DDX3Y exp (TPM)") +
  geom_text( x = 10, y = -0.5, label = "R = 0.2864524")+
  #annotate("label", x = 2.4, y =-3.8 , label = "p value = 0.0006801768") +
  geom_text(label = "p value = 6.799e-06", x = 10, y =-0.55) 


#p-value < 2.2e-16
cor_result = cor.test(cmb_df$DDX3X, cmb_df$DDX3Y)



#EXPRESSION DIFFERENCE of DDX3Y in chrY
DDX3Y_exp = CCLE_exp[,c(1, which(colnames(CCLE_exp) %in% c("DDX3Y")))]
DDX3Y_exp %>% 
  dplyr::filter(DepMap_ID %in% c(Y_loss_sample,Y_normal)) %>% 
  dplyr::mutate(aneu_com = ifelse(DepMap_ID %in% Y_loss_sample, "chrY_loss", "chrY_intact" )) ->P
#left_join(cell_info_flt, by ="DepMap_ID" ) %>% 
#dplyr::filter(Gender %in% c("male", "female")) ->p
P %>% 
  ggboxplot(x = "aneu_com", y = "DDX3Y",
            color = "aneu_com", palette = "jco", fill = "aneu_com",
            alpha = 0.3,
            add = "jitter",width = 0.5) +
  # geom_jitter(alpha = 0.5)+
  geom_rangeframe()+
  theme_tufte() +
  scale_fill_manual(values = c(colors[1],colors[2]))+
  scale_color_manual(values = c(colors[1],colors[2])) +
  labs(x = "", y = "DDX3Y expression")


p_value = t.test(P[P$aneu_com == "chrY_loss",]$DDX3Y, P[P$aneu_com == "chrY_intact",]$DDX3Y)
p_value$p.value

