# HARMONIZE_MATRICES_RSQTL.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 1, 2019

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("tidyverse"); load_package("data.table")

handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # load in matrices
  vaf <<- data.frame(fread(arg_df$value[arg_df$flag == "-r"]))
  spl <<- data.frame(fread(arg_df$value[arg_df$flag == "-s"]))
  cov <<- data.frame(fread(arg_df$value[arg_df$flag == "-c"]))
  
}

# take arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

# get the samples present in each
v.names <- colnames(vaf)[-1]
s.names <- colnames(spl)[-1]
c.names <- colnames(cov)[-1]

# get the samples in all 3 matrices
conc <- intersect(intersect(v.names, s.names), c.names)

# select only the concordant samples from all 3 matrices
vaf.clean <- vaf[,c("SNP", conc)]
spl.clean <- spl[,c("intron", conc)]
cov.clean <- cov[,c("id", conc)]

# write the fixed outputs to a file
write.table(vaf.clean, paste0(gsub(".txt", "", arg_df$value[arg_df$flag == "-r"]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
write.table(spl.clean, paste0(gsub(".txt", "", arg_df$value[arg_df$flag == "-s"]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
write.table(cov.clean, paste0(gsub(".txt", "", arg_df$value[arg_df$flag == "-c"]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
