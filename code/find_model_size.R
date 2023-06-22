## ------------------------------------------------------------------------
## Title: Find optimal model size
## Purpose: Build linear models of various size to find the number of TFs that balances predicitve performance with parsimony
## Author: Roshan Darji
## Date Created: 2021-05-18
## ------------------------------------------------------------------------
## Notes:
##  This script takes 10 terminal arguments (I'm sorry):
##      1. input_data_path - Path to the data
##      2. expr_col - Name of the column with the expression information
##      3. binding_start_col - Name of the first column with TF binding information. All TF binding information needs to be in a contiguous block of columns.
##      4. log_transform - Do you want to log transform the expression information? TRUE or FALSE
##      5. pseudocount - Numeric that is added to the expression column prior to log transformation. Should be set to 0 if log_transform is FALSE.
##      6. num_tf_start - Starting value of the sequence for the number of TFs to include in the linear model
##      7. num_tf_end - Starting value of the sequence for the number of TFs to include in the linear model
##      8. num_tf_by - Incriment of the sequence
##      9. models_output_path - Path to output for model data
##      10. plot_output_path - Path to outupt for plots
##
##  Copy and paste the following, replacing the variables with your own, into the terminal to run this script:
##      Rscript /path/to/find_model_size.R \
##          input_data_path \
##          expr_col \
##          binding_start_col \
##          log_transform \
##          pseudocount \
##          num_tf_start \
##          num_tf_end \
##          num_tf_by \
##          models_output_path \
##          plot_output_path
##
## ------------------------------------------------------------------------



# Packages and necessary funcitons --------------------------------------------------------------------------------
library(tidyverse)
source("/gpfs/gpfs1/home/bmoyers/Scripts/ENCODE_500Plus/Roshan_Cobinding/binding-expression-modeling/rscripts/binding_expr_functions.R")



# Get terminal arguments ------------------------------------------------------------------------------------------
terminal_args <- commandArgs(trailingOnly = TRUE)
# Set up general model parameters
input_data_path <- terminal_args[1] %>% as.character()
expr_col <- terminal_args[2] %>% as.character()
binding_start_col <- terminal_args[3] %>% as.character()
log_transform <- terminal_args[4] %>% as.logical()
pseudocount <- terminal_args[5] %>% as.numeric()

# Set up vector of number of transcription factors
num_tf_start <- terminal_args[6] %>% as.numeric()
num_tf_end <- terminal_args[7] %>% as.numeric()
num_tf_by <- terminal_args[8] %>% as.numeric()

num_tf <- seq.int(from = num_tf_start, to = num_tf_end, by = num_tf_by)

models_output_path <- terminal_args[9]
plot_output_path <- terminal_args[10]

# Build models ----------------------------------------------------------------------------------------------------
models <- map(
    num_tf,
    ~ build_binding_expr_models(input_data_path = input_data_path,
                                expr_col = expr_col,
                                binding_start_col = binding_start_col,
                                log_transform = log_transform,
                                pseudocount = pseudocount,
                                type = "mlr",
                                num_tf = .x,
                                num_cv = 10,
                                num_reps = 100) %>%
        # Only extract the model data, not the summary data
        .$models %>%
        # Each term from the same model will have the same model stats, so we'll filter for the "(Intercept)" term,
        # which is included in every model.
        # This way we only get the model stats once instead of multiple times.
        filter(term == "(Intercept)") %>%
        # Tag on the number of TFs used in the model
        mutate(num_tf = .x)
)


# Save model data -------------------------------------------------------------------------------------------------
saveRDS(models,
        file = models_output_path)


# Plotting --------------------------------------------------------------------------------------------------------
# Transform data to long format to prepare for plotting
models <- models %>%
    bind_rows() %>%
    mutate(num_tf = as.factor(num_tf)) %>%
    pivot_longer(cols = pearson_cor:aic, names_to = "model_statistic")


png(filename = plot_output_path,
    width = 1800, height = 1000, units = "px", res = 150, type="cairo")

models %>%
    ggplot(aes(num_tf, value)) +
    geom_boxplot() +
    facet_wrap(~ statistic, scales = "free", ncol = 1) +
    labs(title = "Subsampled linear model performance",
         x = "Number of predictors",
         y = "Value")

dev.off()
