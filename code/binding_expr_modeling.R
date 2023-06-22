## ------------------------------------------------------------------------
## Title: Binding-Expression Modeling
## Purpose: Build models of Gene Expression ~ TF Binding
## Author: Roshan Darji
## Date Created: 2021-04-07
## ------------------------------------------------------------------------
## Notes:
##   This file takes in multiple terminal arguments:
##   1. input_data_path_arg -- Path to the input data being modeled
##   2. expr_col_arg -- Name of the column with the gene expression information
##   3. binding_start_col -- Name of the first column with TF binding data; all TF binding data needs to be in contiguous columns
##   4. pseudocount_arg -- Sets the number of pseudocounts to be added before log-transforming the gene expression values
##   5. type_arg -- Sets the type of regression; "mlr" for multivariate and "olr" for univariate
##   6. num_tf_arg -- Number of TFs to be used in each model
##   7. seed_arg -- Set the random seed to use
##   8. save_path -- Path to where you want to save the data; needs to end in ".rds"
##
##   If you want more customization over the model building process, use `source("binding-expression-modeling/rscripts/binding_expr_functions.R")`
##   and it's function: `build_binding_expr_models()`. You'll have finer control over how and what gets modeled that way.
##
##   This script is mainly just to simply recreate the methods I took to build the binding-expression models.
##
##   To quickly use this script in the terminal, copy/paste the following command and replace the args with the acutal arguments you want to use:
##      Rscript binding_expr_modeling.R \
##          input_data_path_arg \
##          expr_col_arg \
##          binding_start_col_arg \
##          pseudocount_arg \
##          type_arg \
##          num_tf_arg \
##          seed_arg \
##          save_path
##
##   Do note that you'll probably need a lot of RAM.
##
## ------------------------------------------------------------------------



# Packages --------------------------------------------------------------------------------------------------------
library(tidyverse)
source("/gpfs/gpfs1/home/bmoyers/Scripts/ENCODE_500Plus/Roshan_Cobinding/binding-expression-modeling/rscripts/binding_expr_functions.R")



# Paths to data ---------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE) # Get terminal arguments
input_data_path_arg <- args[1] # Set input data path
expr_col_arg <- args[2] # Set gene expression colname
binding_start_col_arg <- args[3] # Set the first column with TF binding info
pseudocount_arg <- as.numeric(args[4]) # Sets how many pseudocounts should be added when log-transforming
type_arg <- args[5] # Sets the type of regression (multivariate or univariate)
num_tf_arg <- as.numeric(args[6]) # Sets how many TFs are used to model
seed_arg <- as.numeric(args[7])
save_path <- args[8] # Save path needs to end with .rds

# Build models ----------------------------------------------------------------------------------------------------
set.seed(seed_arg)
binding_expr_models <- build_binding_expr_models(
    input_data_path = input_data_path_arg, expr_col = expr_col_arg, binding_start_col = binding_start_col_arg,
    pseudocount = pseudocount_arg,
    type = type_arg, num_tf = num_tf_arg
)


# Save data -------------------------------------------------------------------------------------------------------
saveRDS(binding_expr_models,
        file = save_path)
