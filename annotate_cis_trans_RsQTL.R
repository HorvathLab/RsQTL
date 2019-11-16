# ANNOTATE_CIS_TRANS_RSQTL.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 1, 2019

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("tidyverse"); load_package("data.table"); load_package("GenomicRanges")

handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # load in a RsQTL result file
  rsqtls <<- fread(arg_df$value[arg_df$flag == "-r"])

  # load in the gene locations file
  gene_locs <<- fread(arg_df$value[arg_df$flag == "-g"])
  gene_locs_gr <<- GRanges(gene_locs)
  
  # specify output prefix
  output_prefix <<- arg_df$value[arg_df$flag == "-o"]
  
}

# read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

# annotate which gene the SNP resides in
# classify RsQTLs in which the two members of the pair are in the same gene as cis
# classify all others as trans
res <- rsqtls %>% mutate(new_SNP = SNP, new_gene = gene) %>%
  separate(new_SNP, into = c("chrom", "start", "ref", "alt")) %>%
  mutate(end = start)

res_gr <- GRanges(res)

overlaps <- findOverlaps(res_gr, gene_locs_gr, select = "last", type = "within")
genes_snp <- gene_locs$ensembl_gene[overlaps]

res <- rsqtls %>% mutate(new_SNP = SNP, new_gene = gene) %>%
  separate(new_gene, into = c("chrom", "start", "end")) %>%
  mutate(chrom = gsub("chr", "", chrom))
res_gr <- GRanges(res)
overlaps <- findOverlaps(res_gr, gene_locs_gr, select = "last", type = "within")
genes_intron <- gene_locs$ensembl_gene[overlaps]

res <- rsqtls
res$gene_snp <- genes_snp
res$gene_intron <- genes_intron

res <- res %>% mutate(class = ifelse(gene_snp == gene_intron, "cis", "trans"))
  
# write the results to a file
write.table(res, paste0(output_prefix, "_cistrans_ann.txt"), quote = F, row.names = F, sep = '\t')