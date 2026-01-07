#' Extract Enrichment Terms for Multiple Conditions and Categories
#'
#' This function extracts enrichment terms for multiple conditions and database categories.
#'
#' @param input A data frame containing gene data.
#' @param condition_col A character vector of the column name(s) of the conditions to analyze.
#' @param categories A character vector of database categories for enrichment analysis.
#' @param background A background list of genes to use with enrichR. See enrichr documentation for details.
#' @param split_by_reg Logical. If `TRUE`, splits genes by up- and down-regulation.
#' @param logFC_threshold The threshold for determining up/down-regulation. Defaults to `0`.
#' @param pval_threshold Numeric; adjusted p-value cutoff. Default is 0.05.
#' @param n Integer; number of top terms to return. Default is 10.
#' @param min_genes_threshold Integer; minimum genes to overlap between terms. 
#' @param enrichment_method Which method to use for enrichment, "enrichr" (requires internet connection), or "gmt" (offline, requires gmt file input. 
#' @param gmt_files Required if enrichment_method = "gmt". Path to one or more local gmt file(s). For gmt data format see https://docs.gsea-msigdb.org/#GSEA/Data_Formats/. 
#' @param gmt_min_overlap Minimum overlap between geneset and dataset
#' @param gmt_min_set_size Minimum size of geneset
#' @param gmt_max_set_size Maximum size of geneset
#' @param verbose Logical. If TRUE, print informative messages about the progress. 
#' @param show_progress Logical. If TRUE, display a progress bar during computation.
#' @return A combined data frame of enrichment terms for all specified conditions and categories.
#' @examples
#' # Example usage
#' input <- data.frame(genes = c("Gene1", "Gene2", "Gene3"))
#' result <- extractMultipleTerms(
#'   input, 
#'   condition_col = c("Condition1", "Condition2"), 
#'   categories = c("GO_Biological_Process", "KEGG"), 
#'   top = TRUE
#' )
#' @importFrom dplyr %>% mutate
#' @importFrom tidyr unite
#' @importFrom factoextra fviz_dist
#' @export
extractMultipleTerms <- function(input, condition_col, categories = NULL, background = NULL, split_by_reg = FALSE,
                                 logFC_threshold = 0, pval_threshold = 0.05, n = 10, min_genes_threshold = 3,
                                 enrichment_method = c("enrichr", "gmt"),
                                 gmt_files = NULL,
                                 gmt_min_overlap = 1,
                                 gmt_min_set_size = 5,
                                 gmt_max_set_size = 5000,
                                 verbose = FALSE,
                                 show_progress = interactive()) {
  
  enrichment_method <- match.arg(enrichment_method)
  
  # Check if input is empty
  if (nrow(input) == 0) {
    stop("Error: The input data frame is empty.")
  }
  
  # Check if input contains "genes"
  if ("genes" %in% colnames(input) == FALSE) {
    gene_col <- grep("gene", colnames(input), ignore.case = TRUE, value = TRUE)
    if (length(gene_col) != 1) {
      stop('Error: The input data must contain a gene name column called "genes".')
    }
    colnames(input)[colnames(input) == gene_col] <- "genes"
  }
  
  # Check if input contains "log2fold"
  if (split_by_reg == TRUE && !"log2fold" %in% colnames(input)) {
    stop('Error: The input data must contain column called "log2fold", if split_by_reg is set as true')
  }
  
  # Validate condition_col
  if (!is.character(condition_col) || any(!condition_col %in% colnames(input))) {
    stop("Error: condition_col must be a character vector with valid column names in the input data frame.")
  }
  
  # Combine condition columns into a single composite condition
  input$Composite_Condition <- apply(input[, condition_col, drop = FALSE], 1, paste, collapse = "_")
  
  condition_values <- unique(input$Composite_Condition)
  if (length(condition_values) == 0) {
    stop("Error: No unique values found for the specified condition columns.")
  }
  
  # ---- Categories handling ----
  if (enrichment_method == "enrichr") {
    if (!is.character(categories) || length(categories) == 0) {
      stop("Error: categories must be a valid non-empty character vector for enrichment_method='enrichr'.")
    }
    categories_clean <- categories
    gmt_files_map <- NULL
    
  } else {
    # GMT mode: categories derived from gmt_files (labels) unless provided
    if (is.null(gmt_files) || length(gmt_files) == 0) {
      stop("Error: enrichment_method='gmt' requires 'gmt_files'.")
    }
    
    gmt_files <- as.character(gmt_files)
    
    # Build clean labels for GMT files (strip .gmt from names or basenames)
    if (!is.null(names(gmt_files)) && any(nzchar(names(gmt_files)))) {
      gmt_labels <- sub("\\.gmt$", "", names(gmt_files), ignore.case = TRUE)
    } else {
      gmt_labels <- sub("\\.gmt$", "", basename(gmt_files), ignore.case = TRUE)
    }
    
    # If categories not supplied, use these labels
    if (is.null(categories) || length(categories) == 0) {
      categories_clean <- gmt_labels
    } else {
      categories_clean <- sub("\\.gmt$", "", as.character(categories), ignore.case = TRUE)
    }
    
    # Create a mapping from clean label -> file path
    # (if duplicate labels exist, keep the first and warn)
    if (any(duplicated(gmt_labels))) {
      dup <- unique(gmt_labels[duplicated(gmt_labels)])
      warning("Duplicate GMT labels after stripping '.gmt': ", paste(dup, collapse = ", "),
              ". Using the first occurrence for each.")
      keep_idx <- !duplicated(gmt_labels)
      gmt_files <- gmt_files[keep_idx]
      gmt_labels <- gmt_labels[keep_idx]
    }
    gmt_files_map <- stats::setNames(gmt_files, gmt_labels)
    
    # Validate categories map
    missing_cats <- setdiff(categories_clean, names(gmt_files_map))
    if (length(missing_cats) > 0) {
      stop("Error: These GMT categories could not be matched to gmt_files: ",
           paste(missing_cats, collapse = ", "),
           ". Available: ", paste(names(gmt_files_map), collapse = ", "))
    }
  }
  
  results <- list()
  results_top_terms <- list()
  
  # Progress bar for outer loops
  total_steps <- length(condition_values) * length(categories_clean)
  pb <- NULL
  step <- 0L
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  if (verbose) {
    message("extractMultipleTerms: method=", enrichment_method,
            " | conditions=", length(condition_values),
            " | categories=", length(categories_clean),
            " | total=", total_steps)
  }
  
  for (cond_value in condition_values) {
    subset_input <- input[input$Composite_Condition == cond_value, ]
    
    for (cat in categories_clean) {
      step <- step + 1L
      if (show_progress) utils::setTxtProgressBar(pb, step)
      if (verbose) message("  -> Enriching condition='", cond_value, "' category='", cat, "'")
      
      # Map to per-category GMT file (GMT mode) or leave NULL (Enrichr mode)
      gmt_for_cat <- NULL
      if (enrichment_method == "gmt") {
        gmt_for_cat <- unname(gmt_files_map[[cat]])
      }
      
      filtered_regular_long <- tryCatch(
        readGeneList(
          listpathSig = subset_input,
          name = cond_value,
          background = background,
          category = cat,                 # cat is already ".gmt"-free
          split_by_reg = split_by_reg,
          logFC_threshold = logFC_threshold,
          pval_threshold = pval_threshold,
          n = n,
          min_genes_threshold = min_genes_threshold,
          enrichment_method = enrichment_method,
          gmt_files = gmt_for_cat,         # single GMT file per category
          gmt_min_overlap = gmt_min_overlap,
          gmt_min_set_size = gmt_min_set_size,
          gmt_max_set_size = gmt_max_set_size,
          verbose = verbose,
          show_progress = show_progress
        ),
        error = function(e) {
          warning(paste("Error during enrichment for categories:", cat,
                        "condition:", cond_value, "-", e$message))
          return(NULL)
        }
      )
      
      if (is.null(filtered_regular_long)) next
      
      filtered_regular <- filtered_regular_long[[1]]
      top_list <- filtered_regular_long[[2]]
      
      if (!is.null(filtered_regular)) {
        results[[paste(cond_value, cat, sep = "_")]] <- filtered_regular
        results_top_terms[[paste(cond_value, cat, sep = "_")]] <- top_list
      }
    }
  }
  
  final_results <- do.call(rbind, results)
  final_top_list <- do.call(rbind, results_top_terms)
  
  if (is.null(final_results) || nrow(final_results) == 0) {
    warning("Warning: No results generated for the given categories and conditions.")
  }
  
  if (!is.null(final_results) && nrow(final_results) > 0) {
    final_results <- final_results %>% dplyr::filter(Term %in% unique(final_top_list$Term))
    rownames(final_results) <- NULL
    return(final_results)
  }
  
  invisible(NULL)
}
