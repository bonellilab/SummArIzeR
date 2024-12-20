#' Read Gene List and Perform Enrichment
#'
#' Reads a gene list, optionally splits by regulation, and performs enrichment analysis.
#'
#' @param listpathSig A data frame containing the gene list with columns `genes` and optionally `log2fold`.
#' @param name The name of the condition (used in results).
#' @param category The enrichment database to use.
#' @param split_by_reg Logical. If `TRUE`, splits genes by up- and down-regulation.
#' @param logFC_threshold The threshold for determining up/down-regulation. Defaults to `0`.
#' @param pval_threshold The adjusted p-value threshold for filtering results.
#' @param n The number of top terms to return. Default is 10.
#' @param min_genes_threshold The minimum number of genes required for a term to be kept in the results. Default is 3. 
#' @return A list with two data frames: `regular` (all terms) and `top` (top terms).
#' @import dplyr tidyr enrichR
#' @export
readGeneList <- function(listpathSig, name, category, split_by_reg = FALSE, logFC_threshold = 0, pval_threshold = 0.05, n = 10, min_genes_threshold = 3) {
  if (nrow(listpathSig) == 0) {
    stop("Error: The input gene list (listpathSig) is empty.")
  }
  
  if (!"genes" %in% colnames(listpathSig)) {
    stop("Error: The input gene list must contain a 'genes' column.")
  }
  
  gl <- listpathSig %>% dplyr::rename("ID" = "genes")
  
  # Split by regulation if required
  if (split_by_reg) {
    if (!"log2fold" %in% colnames(gl)) {
      stop("Error: 'log2fold' column is required for splitting by regulation.")
    }
    gl <- gl %>%
      dplyr::mutate(
        regulation = dplyr::case_when(
          log2fold > logFC_threshold ~ "up-regulated",
          log2fold < -logFC_threshold ~ "down-regulated",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(regulation))
  }
  
  # Perform enrichment
  gl_enr <- tryCatch(
    enrichr(gl$ID, databases = category),
    error = function(e) {
      warning(paste("Warning: Enrichment analysis failed for category:", category, "-", e$message))
      return(NULL)
    }
  )
  
  if (is.null(gl_enr) || is.null(gl_enr[[category]]) || nrow(gl_enr[[category]]) == 0) {
    warning(paste("Warning: No enrichment terms found for category:", category))
    return(NULL)
  }
  
  # Process enriched results
  BP_gl <- gl_enr[[category]] %>%
    dplyr::rename(adj_pval = Adjusted.P.value) %>%
    dplyr::select(Term, Genes, adj_pval) %>%
    dplyr::mutate(
      Genes = gsub("\\;", ",", Genes),
      Category = category,
      ID = Term
    )
  
  BP_gl <- BP_gl %>%
    dplyr::mutate(
      dbs = category,
      condition = name
    )
  
  # Filter significant terms
  regular_terms <- BP_gl %>%
    dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
    dplyr::select(Term, Genes, adj_pval, dbs, condition) %>%
    tidyr::separate_rows(Genes, sep = ",") %>%
    dplyr::mutate(Genes = trimws(Genes))
  # Filter out terms with fewer than the minimum number of genes
  regular_terms <- regular_terms %>%
    dplyr::group_by(Term) %>%
    dplyr::filter(n_distinct(Genes) >= min_genes_threshold) %>%
    dplyr::ungroup()
  # Extract top hits (top n terms based on adjusted p-value)
  top_terms <- BP_gl %>%
    rowwise() %>%
    mutate(Gene_Count = str_count(Genes, ",") + 1) %>%  # Count the number of genes
    dplyr::filter(Gene_Count >= min_genes_threshold) %>%                # Filter rows by count
    ungroup() %>%
    dplyr::select(-Gene_Count) %>%                             # Remove helper column 
    dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
    dplyr::slice_min(adj_pval, n = n) %>%
    dplyr::select(Term)
  # Filter regular terms to include only terms present in top_terms
  filtered_regular <- regular_terms %>%
    dplyr::filter(Term %in% unique(top_terms$Term))
  return(filtered_regular)
}


