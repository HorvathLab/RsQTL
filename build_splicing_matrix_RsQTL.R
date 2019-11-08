# BUILD_SPLICING_MATRIX_RSQTL.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 1, 2019

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("tidyverse"); load_package("data.table"); load_package("factoextra")

handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # counts file
  counts <<- fread(arg_df$value[arg_df$flag == "-c"])
  
  # minimum allowable number of reads per cluster
  read_threshold <<- ifelse(length(arg_df$value[arg_df$flag == "--read-thresh"]) > 0, as.numeric(arg_df$value[arg_df$flag == "--read-thresh"]), 30)
  
  # maximum fraction of samples allowed to be NA or heterozygous for a given intron
  na_thresh <<- ifelse(length(arg_df$value[arg_df$flag == "--na-thresh"]) > 0, as.numeric(arg_df$value[arg_df$flag == "--na-thresh"]), 0.2)
  one_thresh <<- ifelse(length(arg_df$value[arg_df$flag == "--one-thresh"]) > 0, as.numeric(arg_df$value[arg_df$flag == "--one-thresh"]), 0.2)
  zero_thresh <<- ifelse(length(arg_df$value[arg_df$flag == "--zero-thresh"]) > 0, as.numeric(arg_df$value[arg_df$flag == "--zero-thresh"]), 0.2)
  
  # specify whether to remove samples with a large proportion of NAs
  remove_high_nas <<- ifelse(length(arg_df$value[arg_df$flag == "--excl-na"]) > 0, as.logical(arg_df$value[arg_df$flag == "--excl-na"]), F)
  na_percentile_thresh <<- ifelse(length(arg_df$value[arg_df$flag == "--perc-thresh"]) > 0, as.numeric(arg_df$value[arg_df$flag == "--perc-thresh"]), 0.9)
  
  # specify whether to remove clusters with only one intron
  remove_singletons <<- ifelse(length(arg_df$value[arg_df$flag == "--excl-sing"]) > 0, as.logical(arg_df$value[arg_df$flag == "--excl-sing"]), F)
  
  # specify output prefix
  output_prefix <<- arg_df$value[arg_df$flag == "-o"]
  
}

# take arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

# Reformat and sum reads with the same coordinates per intron per sample, keeping only one cluster #
counts <- counts %>% mutate(temp = V1) %>%
  separate(temp, into = c("contig", "coord1", "coord2", "cluster"), sep = ':') %>%
  mutate(strand = substr(cluster, nchar(cluster), nchar(cluster)),
         cluster = gsub("_[+-]", "", cluster),
         contig = paste0("chr", contig)) %>%
  select(V1, contig:strand, everything()) %>%
  gather(sample, reads, -V1:-strand) %>%
  select(V1, contig:strand, everything()) %>%
  group_by(contig, coord1, coord2, sample) %>%
  mutate(reads_sum = sum(reads)) %>%
  filter(row_number() == 1) %>%
  ungroup()

# Remove non-chr contigs
counts <- counts %>% filter(grepl("^chr[1-9XY]", contig))

# Compute the sum of number of reads per cluster per sample; 
# Classify clusters with less than the specified fraction of reads to be NA for all introns of the cluster in the corresponding sample
counts <- counts %>% group_by(sample, cluster) %>%
  mutate(cluster_read_sum = sum(reads),
         reads = ifelse(cluster_read_sum < read_threshold, NA, reads)) %>%
  ungroup()

# For each sample, compute the ratio of the number of reads per intron over the sum of reads per cluster
counts <- counts %>% mutate(read_ratio = reads / cluster_read_sum)

# Remove introns that are NA, 1, or 0 in more than the specified fraction of samples 
counts <- counts %>% group_by(contig, coord1, coord2) %>%
  mutate(
    intron = paste0(contig, ":", coord1, "_", coord2),
    frac_na = length(read_ratio[is.na(read_ratio)]) / n(),
    frac_zero = length(read_ratio[read_ratio == 0]) / n(),
    frac_one = length(read_ratio[read_ratio == 1]) / n()
  ) %>% 
  ungroup() %>%
  filter(
    frac_na <= na_thresh,
    frac_zero <= zero_thresh,
    frac_one <= one_thresh
  ) 

# Remove clusters that have only a single intron
if (remove_singletons) {
  counts <- counts %>% group_by(cluster) %>% 
    mutate(num_introns = length(unique(intron))) %>%
    filter(num_introns != 1) %>%
    ungroup()
}

# Make location file from the coordinates on each row
counts_loc <- counts %>% select(intron, contig, coord1, coord2) %>%
  distinct()

# Convert file back to a matrix
counts <- counts %>% select(intron, sample, read_ratio) %>%
  spread(sample, read_ratio)

# Remove samples above the specified %ile for # of NA values 
if (remove_high_nas) {
  na_counts <- colSums(is.na(counts[,-1]))
  high_percentile <- quantile(na_counts, na_percentile_thresh)
  high_na <- names(na_counts[na_counts > high_percentile])
  counts <- counts[,-which(colnames(counts) %in% high_na)]
}

# Write outputs
counts <- as.matrix(counts)
counts_loc <- as.matrix(counts_loc)

write.table(counts, paste0(output_prefix, '_splicing_matrix.txt'), quote = F, row.names = F, sep = '\t')
write.table(counts_loc, paste0(output_prefix, '_splicing-loc_matrix.txt'), quote = F, row.names = F, sep = '\t')
