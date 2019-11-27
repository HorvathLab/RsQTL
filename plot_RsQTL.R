# PLOT_RSQTL.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 26, 2019

# install missing required packages and load packages
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))


handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # load in the required matrix files
  vaf_matrix <<- fread(arg_df$value[arg_df$flag == "-r"]) %>% gather(sample, vaf, -contains("SNV")) %>% drop_na()
  splicing_matrix <<- fread(arg_df$value[arg_df$flag == "-s"]) %>% gather(sample, ratio, -contains("intron")) %>% drop_na()
  
  # load in a results file
  if (length(arg_df$value[arg_df$flag == "-res"]) > 0) {
    results <<- fread(arg_df$value[arg_df$flag == "-res"]) %>% 
      as.data.frame() %>%
      arrange(FDR) %>%
      mutate(pair = paste0(SNP, "_", gene))
  }

  # specify plotting mode
  plot_mode <<- arg_df$value[arg_df$flag == "-m"]
  
  # if bulk plotting, specify number of top correlations to plot
  n_to_plot <<- ifelse((arg_df$value[arg_df$flag == "-n"]) > 0, arg_df$value[arg_df$flag == "-n"], 200)
  n_to_plot <<- ifelse(n_to_plot > nrow(results), nrow(results), n_to_plot)
  
  # if plotting an individual correlation, specify the intron and SNV to plot
  my_snv <<- arg_df$value[arg_df$flag == "-snv"]
  my_intron <<- arg_df$value[arg_df$flag == "-intron"]
  
  # specify output prefix
  output_prefix <<- arg_df$value[arg_df$flag == "-o"]
  
}

# read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

if(plot_mode == "bulk") {
  ### FOR PLOTTING BULK SNVs
  results <- head(results, n_to_plot)
  vaf_matrix <- vaf_matrix[unlist(vaf_matrix[,1]) %in% results[,1],]
  splicing_matrix <- splicing_matrix[unlist(splicing_matrix[,1]) %in% results[,2],]
  
  plot_df <- left_join(left_join(results, vaf_matrix, by = c("SNP" = "SNV")), splicing_matrix, by = c("gene" = "intron", "sample"))
  p <- ggscatter(plot_df, x = "vaf", y = "ratio",
                 fill = "lightsteelblue", color = "lightsteelblue", shape = 21, size = 1.5, # Points color, shape and size
                 add = "reg.line",  # Add regression line
                 add.params = list(color = "dodgerblue4", fill = alpha("dodgerblue4"), 0.5), # Customize reg. line
                 cor.coef = T, # Add correlation coefficient. see ?stat_cor
                 cor.coef.size = 5,
                 cor.coeff.args = list(method = "spearman", label.sep = "\n")) +
    labs(x = "Variant allele fraction", y = "Proportion of read spanning intron junction")
  
  pdf(paste0("output/RsQTL_plots_top", n_to_plot, ".pdf"), width = ifelse(n_to_plot / 8 >= 8, n_to_plot / 8, 8), height = ifelse(n_to_plot / 8 >= 8, n_to_plot / 8, 8))
  print(facet(p, facet.by = "pair", ncol = 5, scales = "free_y"))
  garbage <- dev.off()
  cat(paste0("Bulk plot saved to saved to output/", "RsQTL_plots_top", n_to_plot, ".pdf\n"))
} else if (plot_mode == "single") {
  ### FOR PLOTTING INDIVIDUAL SNVs
  plot_df <- left_join(vaf_matrix %>% filter(SNV == !!my_snv), 
                       splicing_matrix %>% filter(intron == !!my_intron), 
                       by = 'sample')
  
  p <- ggscatter(plot_df, x = "vaf", y = "ratio",
                 fill = "lightsteelblue", color = "lightsteelblue", shape = 21, size = 1.5, # Points color, shape and size
                 add = "reg.line",  # Add regression line
                 add.params = list(color = "dodgerblue4", fill = alpha("dodgerblue4"), 0.5), # Customize reg. line
                 cor.coef = T, # Add correlation coefficient. see ?stat_cor
                 cor.coef.size = 5,
                 cor.coeff.args = list(method = "spearman", label.sep = "\n")) +
    labs(x = "Variant allele fraction", y = "Proportion of read spanning intron junction", title = paste0(my_snv, "-", my_intron))
  
  pdf(paste0("output/", paste0(gsub(":", "_", paste0("RsQTL_plot_", my_snv, "-", my_intron, ".pdf")))))
  print(p)
  garbage <- dev.off()
  cat(paste0("Single plot saved to output/", gsub(":", "_", paste0("RsQTL_plot_", my_snv, "-", my_intron, ".pdf\n"))))
  
} else stop("ERROR: Invalid mode specified!\n")
