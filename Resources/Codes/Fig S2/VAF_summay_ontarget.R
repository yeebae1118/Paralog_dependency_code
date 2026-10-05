# project: calculate the whole length of PPP2CA, PPP2R1A, TYMS and FBXW11, Gene mutation frequency
# Date: April.18.2025
# Author: Yi

# load library
library(tidyverse)
library(ggthemes)
library(ggpubr)

# data
# load data
# list file in VAF folder
File_list = list.files("VAF_folder2",full.names = T)
sample_name = c()
gene_list = c()
for (i in File_list) {
  
  # read each line of the input
  file_name = i
  info_tpm = str_split(file_name, pattern = "/")[[1]][2]
  
  # split the string and get the sample name and gene name 
  gene_name = str_split(info_tpm,pattern = "_")[[1]][[2]]
  sample_name_tpm = paste(str_split(info_tpm,pattern = "_")[[1]][[1]], 
                          str_split(info_tpm,pattern = "_")[[1]][[2]] ,
                          str_split(info_tpm,pattern = "_")[[1]][[3]],sep = "_"
  )
  
  sample_name = c(sample_name,sample_name_tpm)
  gene_list = c(gene_list, gene_name)
}


# make function and read all the VAF files 
READ_VAF = function(x){
  file_name = File_list[x]
  Sample_ID = sample_name[x]
  gene = gene_list[x]
  lines <- readLines(file_name)
  split_lines <- strsplit(lines, "\t")
  parsed <- lapply(split_lines, function(x) {
    chrom <- x[1]
    pos <- as.integer(x[2])
    ref <- x[3]
    depth <- as.integer(x[4])
    
    # Process base fields (can be A, C, G, T, DEL, etc.)
    base_data <- x[5:length(x)]
    base_df <- do.call(rbind, lapply(base_data, function(b) {
      parts <- unlist(strsplit(b, ":"))
      data.frame(
        base = parts[1],
        count = as.integer(parts[2]),
        avg_mq = as.numeric(parts[3]),
        avg_bq = as.numeric(parts[4]),
        stringsAsFactors = FALSE
      )
    }))
    
    base_df$chrom <- chrom
    base_df$pos <- pos
    base_df$ref <- ref
    base_df$depth <- depth
    base_df$sample <- Sample_ID
    base_df$gene <- gene
    
    return(base_df)
  })
  
  result_df <- bind_rows(parsed)
  return(result_df)
}

VAF_list = lapply(c(1:length(File_list)), READ_VAF)
VAF_list_cmb2 = do.call(rbind, VAF_list)

#write_rds(VAF_list_cmb2, "VAF_list_cmb2.RDS")
VAF_list_cmb2 = readRDS("VAF_list_cmb2.RDS")

### check for the entire gene and check the muatation variants on the ontarget events
# the gene I have is PPP2CA, PPP2R1A, FBXW11 and TYMS
paralog_query_ontarget_list = tibble("chr" = c("chr5", "chr19", "chr5", "chr10"),
                                     "start" = c(134194332, 52190052, 171861549, 657653),
                                     "end" = c(134226073, 52229518, 172006638, 673578),
                                     "gene" = c("PPP2CA", "PPP2R1A", "FBXW11","TYMS"))
paralog_query_ontarget_list = as.data.frame(paralog_query_ontarget_list)



ontarget_paralog_freq = function(x){
  # get sample name
  sample_id = sample_name[x]
  
  # get the data for 
  VAF_list_cmb2 %>% 
    dplyr::filter(sample == sample_id) -> df_tpm
  
  # generate the list
  gene_list = c()
  VAF_list = c()
  
  # get the list of gene 
  for (i in c(1:dim(paralog_query_ontarget_list)[1])) {
    # gene name
    gene_name = paralog_query_ontarget_list[i,]$gene
    
    # start site
    start_site = paralog_query_ontarget_list[i,]$start
    
    # end site
    end_site = paralog_query_ontarget_list[i,]$end
    # get the range of gene locus
    df_tpm %>% 
      dplyr::filter(pos >= start_site & pos <= end_site ) -> df_tpm_flt
    
    mutation_list = c()
    total_read = c()
    # get the position size
    pos_list = unique(df_tpm_flt$pos)
    for (i in pos_list) {
      
      df_tpm_flt %>% 
        dplyr::filter(pos == i) -> df_tpm_flt_site
      
      # get the ref 
      ref = df_tpm_flt_site$ref
      
      df_tpm_flt_site %>% 
        filter(base != ref) -> df_tpm_flt_ref
      
      # get the sum of mutant reads
      sum_mutation = sum(df_tpm_flt_ref$count)
      mutation_list = c(mutation_list, sum_mutation)
      # get the total reads 
      total = df_tpm_flt_ref$depth
      total_read = c(total_read,total)
    }
    
    # calculate the frequency of VAF 
    VAF = sum(mutation_list) / sum(total_read)
    
    # add to the list 
    gene_list = c(gene_list, gene_name)
    VAF_list = c(VAF_list, VAF)
  }
  
  return_df = data.frame(sample_id, gene_list, VAF_list)
  return(return_df)
}

paralog_offtarget_cal = lapply(c(1:length(sample_name)), ontarget_paralog_freq)

paralog_offtarget_cal_cmb = do.call(rbind, paralog_offtarget_cal)

