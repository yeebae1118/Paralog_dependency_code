# project: identify the identical sequence in sgRNA in its paralog partner
# Author: Yi
# Date: Yi

# library
library(Biostrings)
library(tidyverse)

# reverse complementary
reverse_complement <- function(seq) {
  comp <- chartr("ACGTacgt", "TGCAtgca", seq)
  return(stringi::stri_reverse(comp))
}

# generate function for the query
Sequence_identity = function(query_seq, ref_seq){
  # Define query and reference
  query <- query_seq
  reference <- ref_seq
  reference = toupper(reference)
  # Initialize
  longest_match <- ""
  query_length <- nchar(query)
  
  # forword
  # Loop through all possible substrings of the query
  for (start_pos in 1:query_length) {
    for (end_pos in start_pos:query_length) {
      subseq <- substr(query, start_pos, end_pos)
      if (grepl(subseq, reference)) {
        # Update if longer than previous match
        if (nchar(subseq) > nchar(longest_match)) {
          longest_match <- subseq
        }
      }
    }
  }
  
  
  # reverse complement
  ref_rev = reverse_complement(reference)
  for(start_pos in 1:query_length){
    for (end_pos in start_pos:query_length) {
      subseq = substr(query, start_pos,end_pos)
      if(grepl(subseq, ref_rev)){
        if(nchar(subseq) >  nchar(longest_match)){
          longest_match <- subseq
        }
      }
    }
  }
  
  # Output the result
  cat("Longest continuous match:", longest_match, "\n")
  cat("Length:", nchar(longest_match), "\n")
  
  ratio = nchar(longest_match) / query_length
  df_return = data.frame(longest_match, ratio)
  return(df_return)
}

# get the list query seq
query_table = readxl::read_xlsx("~/Documents/Postdoc_SheltzerLab/Postdoc_sheltzer/Publication/Paralog dependency/codes/Paralogs identity/query_guide_list.xlsx")

for (i in 1:dim(query_table)[1]) {
  tpm_df = query_table[i,]
  gene = tpm_df$Gene
  query_guide = tpm_df$sgRNA
  ref_info = tpm_df$`reference seq`
  guide = tpm_df$guides
  df_output = cbind(tpm_df[,c(1,2,4,5)],Sequence_identity(query_guide, ref_info))
  
  if(!exists("query_ratio_output")){
    query_ratio_output = df_output
  }else{
    query_ratio_output = rbind(query_ratio_output, df_output)
  }
}


# there is an outlier and some guides mapped into the intron
# to solve the issue, I would need to look into the gDNA
## Load the library
library(biomaRt)

# Connect to Ensembl
mart <- useMart(biomart="ensembl",
                dataset="hsapiens_gene_ensembl")

# Get the DNA sequence of a gene by gene name (e.g., PPP2CB)
seq_data_PPP2CB <- getSequence(
  id = "PPP2CB", 
  type = "hgnc_symbol", 
  seqType = "gene_exon_intron",  # or "gene_flank", "coding", "3utr", etc.
  mart = mart
)

#   longest_match ratio
#  CTCACCTTCTCGCA   0.7
Sequence_identity("TACAGCTCACCTTCTCGCAG",seq_data_PPP2CB$gene_exon_intron )


# Get the DNA sequence of a gene by gene name (e.g., PPP2R1B)
seq_data_PPP2R1B <- getSequence(
  id = "PPP2R1B", 
  type = "hgnc_symbol", 
  seqType = "gene_exon_intron",  # or "gene_flank", "coding", "3utr", etc.
  mart = mart
)

#   longest_match ratio
#   GACACTCGG     0.45
Sequence_identity("TCACAGCACTGGACACTCGG",seq_data_PPP2R1B$gene_exon_intron )

# change the sequence and ratio to certain position of guides
# PPP2CA_giude4
query_ratio_output[4,]$longest_match = "CTCACCTTCTCGCA"
query_ratio_output[4,]$ratio = 0.7

# PPP2R1A_giude4
query_ratio_output[8,]$longest_match = "GACACTCGG"
query_ratio_output[8,]$ratio = 0.45
