#' Extract Enrichment Terms for Multiple Conditions
#'
#' This function extracts enrichment terms for a given input data frame and a list of database categories.
#'
#' @param input A data frame containing gene data.
#' @param condition_col A character vector of the column name of the conditions to analyze.
#' @param category A character vector of database names for enrichment analysis.
#' @param split_by_reg Logical. If `TRUE`, splits genes by up- and down-regulation.
#' @param logFC_threshold The threshold for determining up/down-regulation. Defaults to `0`.
#' @param pval_threshold Numeric; adjusted p-value cutoff. Default is 0.05. 
#' @param n Integer; number of top terms to return. Default is 10.
#'
#' @return A data frame of enrichment terms across all conditions and databases.
#' @examples
#' # Example usage
#' input <- data.frame(genes = c("Gene1", "Gene2", "Gene3"))
#' result <- extract_Terms(input, condition = "Condition1", datab = c("GO_Biological_Process"), top = TRUE)
#' @importFrom dplyr %>% mutate
#' @export
extractTerms <- function(input, condition_col, category, split_by_reg = FALSE, logFC_threshold = 0, 
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
  
  
  # Check if input contains "log2fold"
  if (split_by_reg == T & "log2fold" %in% colnames(input) == F ){
    stop('Error: The input data must contain column called "log2fold", if split_by_reg is set as true')
  }
 
  # Validate condition_col
  if (!is.character(condition_col) || !(condition_col %in% colnames(input))) {
    stop("Error: Condition column must be a valid column name in the input data frame.")
  }
  
  # Check that category is a non-empty character vector
  if (!is.character(category) || length(category) == 0) {
    stop("Error: Category must be a valid non-empty character vector.")
  }
  
  # Extract unique condition values
  condition_values <- unique(input[[condition_col]])
  if (length(condition_values) == 0) {
    stop("Error: No unique values found for the specified condition column.")
  }
  
  results <- list()
  
  # Iterate through condition values
  for (cond_value in condition_values) {
    # Subset input for the current condition
    subset_input <- input[input[[condition_col]] == cond_value, ]
    
    # Iterate through each category in the vector
    for (cat in category) {
      filtered_regular <- tryCatch(
        readGeneList(
          listpathSig = subset_input,
          name =cond_value, 
          category = cat,
          split_by_reg = split_by_reg,
          logFC_threshold = logFC_threshold,
          pval_threshold = pval_threshold,
          n = n,
          min_genes_threshold = min_genes_threshold
        ),
        error = function(e) {
          warning(paste("Error during enrichment for category:", cat, 
                        "condition:", cond_value, "-", e$message))
          return(NULL)
        }
      )
      
      # Store results
      if (!is.null(filtered_regular)) {
        results[[paste(cond_value, cat, sep = "_")]] <- filtered_regular
      }
    }
  }
  
  # Combine results across all conditions and categories
  final_results <- do.call(rbind, results)
  
  # Warning if no results are generated
  if (is.null(final_results) || nrow(final_results) == 0) {
    warning("Warning: No results generated for the given categories and conditions.")
  }
  
  return(final_results)
}


