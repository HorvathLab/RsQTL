# BUILD_PCA_COVARIATE_MATRIX_RSQTL.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 26, 2019

# load packages
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(factoextra))

handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # splicing and read count matrix files
  rc <<- data.frame(fread(arg_df$value[arg_df$flag == "-r"]))
  splicing <<- data.frame(fread(arg_df$value[arg_df$flag == "-s"]))
  additional_covs <<- data.frame(ifelse(length(arg_df$value[arg_df$flag == "-c"]) > 0, fread(arg_df$value[arg_df$flag == "-c"]), data.frame()))
  
  # specify how many PCs to retain
  n_pcs <<- ifelse(length(arg_df$value[arg_df$flag == "-n"]) > 0, as.numeric(arg_df$value[arg_df$flag == "-n"]), 10)
  
  # specify output prefix
  output_prefix <<- arg_df$value[arg_df$flag == "-o"]
  
}

# take arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

# compute pcs for splicing data
pcs_spl <- prcomp(data = splicing[,-1], ~., na.action = na.omit)
pca_spl_plot <- fviz_eig(pcs_spl, addlabels = TRUE)

# output pcs plot for splicing data
pdf(paste0("output/", output_prefix, "_splicing_pca_plot.pdf"))
print(pca_spl_plot)
garbage <- dev.off()
cat(paste0("Splicing PCA plot saved to output/", output_prefix, "_splicing_pca_plot.pdf\n"))

# compute pcs for read count data
pcs_rc <- prcomp(data = rc[,-1], ~., na.action = na.omit)
pca_rc_plot <- fviz_eig(pcs_rc, addlabels = TRUE)

# output pcs plot for readcount data
pdf(paste0("output/", output_prefix, "_readcount_pca_plot.pdf"))
print(pca_rc_plot)
garbage <- dev.off()
cat(paste0("Read-count PCA plot saved to output/", output_prefix, "_readcount_pca_plot.pdf\n"))

# select specified number of pcs (if possible)
n_pcs_act <- min(n_pcs, ncol(pcs_spl$rotation), ncol(pcs_spl$rotation))
if(n_pcs_act != n_pcs) cat(paste0("Too few observations to compute ", n_pcs, " PCs. Using the maximum possible number (", n_pcs_act, ") instead.\n"))
pcs_spl <- pcs_spl$rotation[,1:n_pcs_act]
pcs_rc <- pcs_rc$rotation[,1:n_pcs_act]

# combine into one file
pcs_spl <- data.frame(t(pcs_spl))
pcs_spl$id <- paste0("SPLICING_", rownames(pcs_spl))

pcs_rc <- data.frame(t(pcs_rc))
pcs_rc$id <- paste0("RC_", rownames(pcs_rc))

covs_all <- bind_rows(pcs_spl, pcs_rc, additional_covs) %>% select(id, everything())

# output covariate matrix
write.table(covs_all, paste0("output/", output_prefix, '_covariate_pca_matrix.txt'), quote = F, row.names = F, sep = '\t')
cat(paste0("Covariate matrix saved to output/", output_prefix, "_covariate_pca_matrix.txt\n"))
