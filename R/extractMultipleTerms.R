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
#' 
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
extractMultipleTerms <- function(input, condition_col, categories, background = NULL, split_by_reg = FALSE, logFC_threshold = 0, 
                         pval_threshold = 0.05, n = 10, min_genes_threshold = 3) {

  
  # Check if input is empty
  if (nrow(input) == 0) {
    stop("Error: The input data frame is empty.")
  }
  # Check if input contains "genes"
  if ("genes" %in% colnames(input) == F) {
   gene_col<- grep("gene", colnames(input), ignore.case = T, value = TRUE)
   if(length(gene_col) != 1){
     stop('Error: The input data must contain a gene name column called "genes".')}
   if(length(gene_col) == 1){
     colnames(input)[colnames(input) == gene_col] <- "genes"
     
   }
  }
  
  
  # Check if input is contains "log2fold"
  if (split_by_reg == T & "log2fold" %in% colnames(input) == F ){
    stop('Error: The input data must contain column called "log2fold", if split_by_reg is set as true')
  }
 
  # Validate condition_col as a character vector of valid column names
  if (!is.character(condition_col) || any(!condition_col %in% colnames(input))) {
    stop("Error: condition_col must be a character vector with valid column names in the input data frame.")
  }
  
  # Check that categories is a non-empty character vector
  if (!is.character(categories) || length(categories) == 0) {
    stop("Error: categories must be a valid non-empty character vector.")
  }
  
  # Combine condition columns into a single composite condition
  input$Composite_Condition <- apply(input[, condition_col, drop = FALSE], 1, paste, collapse = "_")
  
  # Extract unique composite condition values
  condition_values <- unique(input$Composite_Condition)
  if (length(condition_values) == 0) {
    stop("Error: No unique values found for the specified condition columns.")
  }
  
  results <- list()
  results_top_terms<-list()
  # Iterate through composite condition values
  for (cond_value in condition_values) {
    # Subset input for the current composite condition
    subset_input <- input[input$Composite_Condition == cond_value, ]
    
    # Iterate through each categories in the vector
    for (cat in categories) {
      filtered_regular_long <- tryCatch(
        readGeneList(
          listpathSig = subset_input,
          name = cond_value, 
	  background = background,
          category = cat,
          split_by_reg = split_by_reg,
          logFC_threshold = logFC_threshold,
          pval_threshold = pval_threshold,
          n = n,
          min_genes_threshold = min_genes_threshold
        ),
        error = function(e) {
          warning(paste("Error during enrichment for categories:", cat, 
                        "condition:", cond_value, "-", e$message))
          return(NULL)
        }
      )
      filtered_regular<-filtered_regular_long[[1]]
      top_list<-filtered_regular_long[[2]]
      # Store results
      if (!is.null(filtered_regular)) {
        results[[paste(cond_value, cat, sep = "_")]] <- filtered_regular
      }
      if (!is.null(filtered_regular)) {
        results_top_terms[[paste(cond_value, cat, sep = "_")]] <- top_list
      }
      
    }

    
  }
  
  # Combine results across all conditions and categories
  final_results <- do.call(rbind, results)
  # Combine top terms across all conditions and categories
  final_top_list <- do.call(rbind, results_top_terms)
 # Warning if no results are generated
  if (is.null(final_results) || nrow(final_results) == 0) {
    warning("Warning: No results generated for the given categories and conditions.")
  }
  # Filter for top hits only 
  if (nrow(final_results) > 0) {
     final_results<-final_results %>% dplyr::filter(Term %in% unique(final_top_list$Term))
  rownames(final_results)<- NULL
      return(final_results)
  }
}




