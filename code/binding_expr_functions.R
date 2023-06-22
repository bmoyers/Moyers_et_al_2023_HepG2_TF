## ------------------------------------------------------------------------
## Title: Expression ~ Binding Modeling Functions
## Purpose:
## Author: Roshan Darji
## Date Created: 2021-04-06
## ------------------------------------------------------------------------
## Notes:
## This is a group of functions based on the scripts I wrote build and summarize the linear models earlier.
##
## Functionally, all the steps from the earlier functions are nearly identical. I removed some of the labeling because
## that data wouldn't always be present/necessary, and can be manually added in later (like data type and genomic region).
## Outside of that the core model building is copy/pasted from my old scripts and small tweaks were made to get it to fit
## properly in these new functions.
##
## For fitting the models, instead of running the modeling through a `for` loop (which are super slow in R), I took the
## contents of the loop and turned it into a function. Models are then created in a vectorized fashion using the `map`
## family of functions from the `purrr` package. This should lead to faster model building.
##
##
##
## ------------------------------------------------------------------------


# Packages required -------------------------------------------------------
# I really lean on the tidyverse suite of packages to manipulate data, so certain packages needed to be available to R.
# The following packages need to be present in `installed.packages()`
# - vroom
# - janitor
# - tidyverse suite (https://www.tidyverse.org/packages/)
# - modelr
# - yardstick
# - GeneBook

# The magrittr and rlang packages needs to be loaded to use the pipe operator (%>%) and quasiquotation in functions
library(magrittr) # If this ever gets added to a package, drop this library() statment and use usethis::use_pipe() in your package to add the pipe functionality to the package.
library(rlang) # Allows for quasiquotation
library(progress) # Progress bars


# Functions ---------------------------------------------------------------
#' Preprocess TF Binding Data
#'
#' @description
#' Prepares transcription factor binding and gene expression data for modeling.
#' Both the TF binding information and gene expression information needs to be within the same dataset, and columns
#' containing TF binding information need to be in a contiguous block (i.e. they can be referenced like TF_1:TF_n).
#'
#'
#' @param input_data Input data frame with TF binding/expression data.
#' @param input_data_path Path to the input TF binding/expression data from the current working directory.
#'     Using absolute paths (starting from root) is recommended to prevent issues.
#' @param expr_col A character string giving the name of column with the gene expression data
#' @param binding_start_col A character string giving the name of the first column with transcription factor binding data
#' @param binding_end_col A character string giving the name of the last column with transcription factor binding data.
#'     If a character string is not provided, then it will assume that the last column with binding data is the last
#'     column in the data frame.
#' @param log_transform A logical value that indicates if `expr_col` needs to be log-transformed prior to use in the
#'     model. Defaults to TRUE.
#' @param pseudocount A numeric value that is added to the gene expression column prior to log-transformation.
#'     Defaults to 0.
#'
#' @return A data frame with the bare minimum data necessary for modeling
#'
preprocess_binding_expr <- function(input_data = NULL, input_data_path = NULL, expr_col, binding_start_col,
                                    binding_end_col = NULL, log_transform = TRUE, pseudocount = 0) {
    # Check if you have the right packages for this
    if (all(!(c("vroom", "janitor", "dplyr", "tidyr", "tidyselect") %in% installed.packages()[, 1]))) {
        stop("You don't have the right packages available. Please make sure that vroom, janitor, dplyr, tidyr, and tidyselect are installed and available to use.")
    }

    # Load input data
    if (!is.null(input_data) && is.null(input_data_path)) {
        input_data <- input_data %>%
            janitor::clean_names() # Cleans up column names. Resulting names are unique and consist of only the "_" character, numbers, and lowercase letters.
    } else if (is.null(input_data) && !is.null(input_data_path)){
        message("Loading input data with vroom \n")
        #input_data <- vroom::vroom(input_data_path) %>% # `vroom` package is super fast at loading data
        input_data <- read.delim(input_data_path, header=T, sep="\t", stringsAsFactors=F) %>% # `vroom` package is super fast at loading data
            janitor::clean_names() # Cleans up column names. Resulting names are unique and consist of only the "_" character, numbers, and lowercase letters.
    } else if (is.null(input_data) && is.null(input_data_path)) {
        stop("Please provide either the input data as a data frame (or equivalent; i.e. tibble) or a path to the input data.")
    } else if (!is.null(input_data) && !is.null(input_data_path)) {
        stop("Please ensure that either input_data or input_data_path is specified, not both.")
    }

    # Reformat expr_col, binding_start_col, and binding_end_col to match the same format as in janitor::clean_names() and covert to symbols
    expr_col <- janitor::make_clean_names(expr_col) %>% sym()
    binding_start_col <- janitor::make_clean_names(binding_start_col) %>% sym()

    if (!is.null(binding_end_col)) {
        binding_end_col <- janitor::make_clean_names(binding_end_col) %>% sym()
    }

    # Remove any NAs from the expr_col and log-transform it (if option is set to TRUE)
    if (log_transform) {
        input_data <- input_data %>%
            dplyr::filter(!is.na(!!expr_col)) %>%
            dplyr::mutate(expr = log(!!expr_col + pseudocount)) %>% # Calculate log(expr_col); pseudocount added to avoid -Inf
            dplyr::relocate(expr, .after = !!expr_col) # Moves the new column to be after expr_col
    } else {
        input_data <- input_data %>%
            dplyr::filter(!is.na(!!expr_col)) %>%
            dplyr::mutate(expr = !!expr_col) %>%
            dplyr::relocate(expr, .after = !!expr_col) # Moves the new column to be after expr_col
    }

    # Tranform the data to long format
    if (is.null(binding_end_col)) { # If binding_end_col is not specified, then use the last column in the data frame
        preprocessed_data <- input_data %>%
            tidyr::pivot_longer(cols = !!binding_start_col:tidyselect::last_col(),
                                names_to = "tf", names_repair = "minimal")
    } else {
        preprocessed_data <- input_data %>%
            tidyr::pivot_longer(cols = !!binding_start_col:!!binding_end_col,
                                names_to = "tf", names_repair = "minimal")
    }

    return(preprocessed_data)
}

#' Detailed progress tracker
#'
#' @description
#' A simple helper function to update the progress bar. Helps to make some of the code later on easier to read.
#'
#' @param use A logical indicating whether or not detailed_progress should be used.
update_detailed_progress <- function(use) {
    # Add a tick to the progress bar if detailed_progress option is used
    if (use) {
        pb$tick()
    }
}

#' Fit a linear regression of gene expression using transcription factor binding information as predictors
#'
#' @description
#' This function takes in the preprocessed data from `preprocess_binding_expr()` and builds a linear regression.
#' Can be run either as a univariate, ordinary linear regression (using `type = "olr"`) or as a multivariate, multiple
#' linear regression (using `type = "mlr"`).
#' For OLR models, a single transcription factor's binding information is used as the predictor for gene expression.
#' For MLR models, a character vector of transcription factors is provided and the model samples this vector for the
#' number of transcription factors specified (set with `num_tf`).
#' If you want to use all transcription factors listed in the input, set the number of transcription factors to be used
#' equal to the length of the input vector.
#'
#' The model creates a number of Monte Carlo cross-validation sets to be used for training and testing.
#' These sets are generated by randomly assigning rows from the model data to be either in the training set or the
#' testing set.
#' The proportion of samples to be in the testing set is specified by `test_prop`.
#' Combining the training and testing sets generated will result in a copy of the original data set.
#'
#'
#' @param preprocessed_data Data frame generated by `preprocess_binding_expr()`
#' @param input_tf A character string or vector with the transcription factor(s) that are to be sampled from as predictors
#' @param type "mlr" for multiple linear regression or "olr for ordinary linear regression
#' @param detailed_progress A logical indicated whether or not a more detailed progress bar should be used.
#'     As this function does not initialize the progress bar from the `progress` package. This needs to be done first,
#'     under the variable name `pb`. Defaults to `FALSE`.
#' @param num_tf The number of transcription factors to be included in the model. Defaults to 45.
#' @param num_cv The number of Monte Carlo cross-validation sets used to train the model. Defaults to 100.
#' @param test_prop Sets the testing data proportion for the train/test split. Defaults to 0.3.
#'
#' @return A data frame containing model parameters from the trained linear models and performance statistics of the
#'     models on the test data.

fit_binding_expr_lm <- function(preprocessed_data, input_tf,
                                type = c("mlr", "olr"), detailed_progress = FALSE,
                                num_tf = 45, num_cv = 100, test_prop = 0.3) {
    # Check if you have the right packages for this
    if (all(!(c( "dplyr", "tidyr", "broom", "yardstick", "modelr") %in% installed.packages()[, 1]))) {
        stop("You don't have the right packages available. Please make sure that broom, dplyr, modelr, tidyr, and yardstick are installed and available to use.")
    }
    # If detailed_progress is being used, check if a progress bar has been initialized.
    if (detailed_progress) {
        if (!("pb" %in% ls(envir = .GlobalEnv))) {
            stop("A progress bar named `pb` has not been initalized in the global environment.")
        }
    }

    # If type is not set, set type to "mlr"
    type <- type[1]
    # Subset for sampled TFs
    if(type == "olr") {
        # Check if only 1 input TF was specified
        if (length(input_tf) != 1) {
            stop("For OLS modeling, input_tf needs to be a single character string of the target TF.")
        }

        # Add a tick to the progress bar if detailed_progress option is used
        update_detailed_progress(detailed_progress)

        # Preprocess data
        model_data <- preprocessed_data[preprocessed_data$tf == input_tf, ] %>%  # Subset for input_tf
            tidyr::pivot_wider(names_from = tf, values_from = value) %>%  # Convert back to a wider data frame for modeling
            dplyr::select(expr, !!sym(input_tf)) # Selects only the columns with gene expression data and specified input transcription factors

    } else if (type == "mlr") {
        # Check if multiple input TFs were specified
        if (!(length(input_tf) > 1)) {
            stop("For MLS modeling, input_tf needs to be a character vector of the target TFs.")
        }
        # Check if number of input TFs is less than the number of TFs specified to be in the model
        if (length(input_tf) < num_tf) {
            stop(paste("Please make sure that the number of TFs in the input is greater than the number of TFs specified to be in the model.\n",
                       num_tf, "TFs specified to be in the model, but only", length(input_tf), "input TFs found."))
        }

        # Add a tick to the progress bar if detailed_progress option is used
        update_detailed_progress(detailed_progress)

        # Preprocess data
        sampled_tfs <- sample(input_tf, num_tf) # Randomly sample TFs
        model_data <- preprocessed_data[preprocessed_data$tf %in% sampled_tfs, ] %>%  # Subset for sampled_tfs
            tidyr::pivot_wider(names_from = tf, values_from = value) %>%  # Convert back to a wider data frame for modeling
            dplyr::select(expr, !!!syms(sampled_tfs)) # Selects only the columns with gene expression data and sampled input transcription factors
    } else {
        stop('Please make sure type is either "mlr" or "olr".')
    }

    # Model prep and training
    model_data.cv <- modelr::crossv_mc(model_data, n = num_cv, test = test_prop);  # Create train/test sets
    model <- purrr::map(model_data.cv$train, ~ lm(expr ~ ., data = .x)) # Fit models
    update_detailed_progress(detailed_progress)
    # Collect model metrics
    model_tidy <- purrr::map(model, ~ broom::tidy(.x)) # Get model summary stats in a data frame
    model_aug <- purrr::map2(model, model_data.cv$test, ~ broom::augment(.x, newdata = .y)) # Predict on test data
    model_metrics <- purrr::map2( # Calculate model metrics
        model_aug, model,
        ~ {
            dplyr::summarize(
                .data = .x,
                pearson_cor = cor(expr, .fitted, method = "pearson"), # Pearson correlation
                spearman_cor = cor(expr, .fitted, method = "spearman"), # Spearman correlation
                yardstick::rsq_trad(.x, truth = expr, estimate = .fitted) %>% dplyr::transmute(rsq = .estimate), # R-squared
                yardstick::rmse(.x, truth = expr, estimate = .fitted) %>% dplyr::transmute(rmse = .estimate)# RMSE
            ) %>%
                dplyr::mutate(aic = AIC(.y)) # AIC
        }
    )
    update_detailed_progress(detailed_progress)

    # Modeling output
    model_output_loop <- purrr::map2_df(model_tidy, model_metrics, ~ dplyr::bind_cols(.x, .y)); update_detailed_progress(detailed_progress) # Bind together model summary stats and metrics
    # Add in what TF(s) we looked at with this model
    if (type == "olr") { # For OLS models this is the input_tf
        model_output_loop <- model_output_loop %>%
            dplyr::mutate(tf = input_tf)
    } else if (type == "mlr") { # For MLS models this is what TFs we sampled
        model_output_loop <- model_output_loop %>%
            dplyr::mutate(tf = list(sampled_tfs))
    }

    return(model_output_loop)
}


#' Summarize group of linear regression models
#'
#' @description
#' Takes in data from `fit_binding_expr_lm()`, and creates a summary dataset from that information.
#' For each parameter in the linear regression, the mean and standard deviation for all statistics are calculated.
#' Fraction of times a term is significant when included in the model is also calculated by counting the number of times
#' a term is included in a model and comparing it to the number of times it appears after filtering for significance.
#' (NOTE: This will always evaluate to 1 if you do not assign a p.value_filter or p.adj_filter).
#' Gene information is also joined to this summary data using the GeneBook::genecard_description_summary dataset.
#' This will not be done if you do not have the GeneBook package installed and available to use.
#'
#'
#' @param model_info Data frame generated by `fit_binding_expr_lm()`.
#' @param p.adj_method Sets the method to be used to adjust p-values. A character value from `p.adjust.methods`.
#'     The default is set to "fdr".
#' @param p.value_filter Filters model parameters whose p-values are less than the specified value.
#' @param p.adj_filter Filters model parameters whose adjusted p-values are less than the specified value.
#'
#' @return A data frame of summarized model parameters

summarize_binding_expr_lm <- function(model_info,
                                      p.adj_method = p.adjust.methods,
                                      p.value_filter = NULL, p.adj_filter = NULL) {
    # If p.adj_method is not set, then set it to "fdr"
    if (missing(p.adj_method)) {
        p.adj_method <- "fdr"
    }

    # Add adjusted p-values to model information
    model_info <- model_info %>%
        dplyr::mutate(p.adj = p.adjust(p.value, p.adj_method)) %>%
        dplyr::relocate(p.adj, .after = p.value)
    # This creates the total number of models run with a specific term
    # This does NOT look in the TF column, but should not need to in order to get an accurate total.
    merged_summary <- model_info %>%
        dplyr::group_by(term) %>%
        dplyr::summarize(n = dplyr::n()) %>%
        dplyr::mutate(id = glue::glue("{term}"))

    # Filter data based on a p.value cut-off
    if (!is.null(p.value_filter)) {
        model_info <- model_info %>% dplyr::filter(p.value < p.value_filter)
    }
    # Filter data based on p.adj cut-off
    if(!is.null(p.adj_filter)) {
        model_info <- model_info %>% dplyr::filter(p.adj < p.adj_filter)
    }

    # Summarize data using filtered data
    models_summarized <- model_info %>%
        dplyr::group_by(term) %>% # Group the data by data_type
        dplyr::summarize(n = dplyr::n(), # This gets the count for each unique combination of the grouped variables above
                  mean_est = mean(estimate),
                  sd_est = sd(estimate), # This is a really naive standard deviation estimate
                  mean_pearson_cor = mean(pearson_cor),
                  sd_pearson_cor = sd(pearson_cor),
                  mean_spearman_cor = mean(spearman_cor),
                  sd_spearman_cor = sd(spearman_cor),
                  mean_rsq = mean(rsq),
                  sd_rsq = sd(rsq),
                  mean_rmse = mean(rmse),
                  sd_rmse = sd(rmse),
                  mean_aic = mean(aic),
                  sd_aic = sd(aic)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(id = glue::glue("{term}"))
    # Using merged_summary from earlier, calculate the number of significant parameters out of the total number of models
    models_summarized <- models_summarized %>%
        dplyr::mutate(frac_sig = n/{merged_summary %>% dplyr::filter(id %in% models_summarized$id) %>% .$n}) %>%
        dplyr::relocate(frac_sig, .after = n) %>%
        dplyr::select(-id) %>%
        dplyr::arrange(-frac_sig, -n)
    # Add in gene descriptions if the GeneBook package is available
    if ("GeneBook" %in% installed.packages()[, 1]) {
        gene_info <- GeneBook::genecard_description_summary %>%
            tibble::as_tibble() %>%
            dplyr::mutate(gene = stringr::str_remove(gene, '["]*')) %>%
            dplyr::mutate(gene = stringr::str_remove(gene, '[\\\\]*"')) %>%
            dplyr::mutate(gene = stringr::str_to_lower(gene)) %>%
            dplyr::select(gene, type, description, summary_genecard)

        models_summarized <- models_summarized %>%
            dplyr::mutate(gene = stringr::str_remove(term, "_flag")) %>%
            dplyr::left_join(gene_info)

    }

    return(models_summarized)
}


#' Build linear regressions of gene expression using transcription factor binding information
#'
#' @description
#' Work-flow for fitting linear regressions of gene expression using transcription factor binding information as predictors.
#' Combines `preprocess_binding_expr()`, `fit_binding_expr_lm()`, and `summarize_binding_expr_lm()` to build large sets
#' of models at one time.
#'
#' @param input_data Input data frame with TF binding/expression data.
#' @param input_data_path Path to the input TF binding/expression data from the current working directory.
#'     Using absolute paths (starting from root) is recommended to prevent issues.
#' @param expr_col A character string giving the name of column with the gene expression data
#' @param binding_start_col A character string giving the name of the first column with transcription factor binding data
#' @param binding_end_col A character string giving the name of the last column with transcription factor binding data.
#'     If a character string is not provided, then it will assume that the last column with binding data is the last
#'     column in the data frame.
#' @param log_transform A logical value that indicates if `expr_col` needs to be log-transformed prior to use in the
#'     model. Defaults to TRUE.
#' @param pseudocount A numeric value that is added to the gene expression column prior to log-transformation.
#'     Defaults to 0.
#' @param input_tf A character string or vector with the transcription factor(s) that are to be sampled from as predictors.
#' @param type "mlr" for multiple linear regression or "olr for ordinary linear regression
#' @param detailed_progress Option to use a more detailed progress bar. The bar will update as parts of the modeling
#'     function finish. Defaults to `TRUE.`
#' @param num_tf The number of transcription factors to be included in the model. Defaults to 45.
#' @param num_cv The number of Monte Carlo cross-validation sets used to train the model. Defaults to 100.
#' @param num_reps The number of samples of transcription factors to be generated (only used in MLR models). Defaults to 500
#' @param test_prop Sets the testing data proportion for the train/test split. Defaults to 0.3.
#' @param p.adj_method Sets the method to be used to adjust p-values. A character value from `p.adjust.methods`.
#'     The default is set to "fdr".
#' @param p.value_filter Filters model parameters whose p-values are less than the specified value.
#' @param p.adj_filter Filters model parameters whose adjusted p-values are less than the specified value. Defaults to 0.05.
#'
#' @return A list of 2 data frames. The first is the raw model parameters and performance statistics for each model.
#'     The second is the summarized model parameters and performance statistics for each transcription factor.

build_binding_expr_models <- function(input_data = NULL, input_data_path = NULL, expr_col, binding_start_col, # Necessary arguments
                                      binding_end_col = NULL, log_transform = TRUE, pseudocount = 0, # Preprocessing arguments
                                      input_tf = NULL, type = c("mlr", "olr"), detailed_progress = TRUE, num_tf = 45, num_cv = 100, test_prop = 0.3, num_reps = 500 ,# Modeling arguments
                                      p.adj_method = p.adj_method, p.value_filter = NULL, p.adj_filter = 0.05) { # Summarization arguments

    # If type is not set, set type to "mlr"
    type <- type[1]
    # If p.adj_method is not set, then set it to "fdr"
    if (missing(p.adj_method)) {
        p.adj_method <- "fdr"
    }

    # Preprocess data
    preproc <- preprocess_binding_expr(input_data, input_data_path, expr_col, binding_start_col, binding_end_col, log_transform, pseudocount)

    # Automatically select transcription factors to be put into the model if not provided with an input or list of input transcription factors.
    if (is.null(input_tf)) { # If the input_tf option is null,
        input_tf <- unique(preproc$tf) # the input is set to be a vector of unique transcription factors present in the original dataset
    }
    # Fit regressions
    message(paste0("\n", paste0(rep("-", options()$width - 2), collapse = ""), "\nFitting linear regressions"))
    if (type == "olr") {
        # Progress bar set up; requires global assignment  (hence the <<- operator)
        pb <<- progress_bar$new(total = length(input_tf) * (detailed_progress * 5), # Will only multiply by 5 if detailed_progress is TRUE. This is a faster way to evaluate `if` statements
                               clear = FALSE,
                               format = "   modeling [:bar] :percent eta: :eta   elapsed: :elapsed")
        pb$tick(0) # Show progress bar immediately
        # Build models
        model_info <- purrr::map_dfr(input_tf, ~ {pb$tick(); fit_binding_expr_lm(preproc, .x, type, detailed_progress, num_tf, num_cv, test_prop)})

    } else if (type == "mlr") {
        # Generate samples of TFs
        sampled_tfs <- replicate(num_reps, sample(input_tf, num_tf), simplify = FALSE)
        # Progress bar set up; requires global assignment  (hence the <<- operator)
        pb <<- progress_bar$new(total = length(sampled_tfs) * (detailed_progress * 5), # Will only multiply by 5 if detailed_progress is TRUE. This is a faster way to evaluate `if` statements
                               clear = FALSE,
                               format = "   [:bar] :percent eta: :eta   elapsed: :elapsed")
        pb$tick(0) # Show progress bar immediately
        # Build models
        model_info <- purrr::map_dfr(sampled_tfs, ~ {pb$tick(); fit_binding_expr_lm(preproc, .x, type, detailed_progress, num_tf, num_cv, test_prop)})
        #model_info <- fit_binding_expr_lm(preproc, .x, type, detailed_progress, num_tf, num_cv, test_prop)
    }

    # Summarize model data
    message(paste0("\n", paste0(rep("-", options()$width - 2), collapse = ""), "\nSummarizing models"))
    model_summary <- summarize_binding_expr_lm(model_info, p.adj_method, p.value_filter = p.value_filter, p.adj_filter = p.adj_filter)

    # Combine data in a list form
    model_list <- list(models = model_info,
                       summary = model_summary)

    message(paste0("\n", paste0(rep("-", options()$width - 2), collapse = ""), "\nDone! Finished @ ", Sys.time()))
    return(model_list)
}
