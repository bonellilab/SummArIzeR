#' Read Gene List and Perform Enrichment
#'
#' Reads a gene list, optionally splits by regulation, and performs enrichment analysis.
#'
#' @param listpathSig A data frame containing the gene list with columns `genes` and optionally `log2fold`.
#' @param name The name of the condition (used in results).
#' @param category The enrichment database to use.
#' @param background A background list of genes to use with enrichR. See enrichr documentation for details. 
#' @param split_by_reg Logical. If `TRUE`, splits genes by up- and down-regulation.
#' @param logFC_threshold The threshold for determining up/down-regulation. Defaults to `0`.
#' @param pval_threshold The adjusted p-value threshold for filtering results.
#' @param n The number of top terms to return. Default is 10.
#' @param min_genes_threshold The minimum number of genes required for a term to be kept in the results. Default is 3. 
#' @return A list with two data frames: `regular` (all terms) and `top` (top terms).
#' @import dplyr tidyr
#' @import stringr
#' @export

readGeneList <- function(listpathSig, name, category, background = NULL, split_by_reg = FALSE, logFC_threshold = 0, pval_threshold = 0.05, n = 10, min_genes_threshold = 3) {
  if (nrow(listpathSig) == 0) {
    stop("Error: The input gene list (listpathSig) is empty.")
  }
  
  if (!"genes" %in% colnames(listpathSig)) {
    stop("Error: The input gene list must contain a 'genes' column.")
  }
  
  gl <- listpathSig %>% dplyr::rename("ID" = "genes")
  
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
    
    results <- list()
    top_term_list<-list()
    
    for (reg in c("up-regulated", "down-regulated")) {
      gl_subset <- gl %>% dplyr::filter(regulation == reg)
      
      gl_enr <- tryCatch(
        enrichr(gl_subset$ID, databases = category, background = background),
        error = function(e) {
          warning(paste("Warning: Enrichment analysis failed for", reg, "genes in category:", category, "-", e$message))
          return(NULL)
        }
      )
      
      if (is.null(gl_enr) || is.null(gl_enr[[category]]) || nrow(gl_enr[[category]]) == 0) {
        warning(paste("Warning: No enrichment terms found for", reg, "genes in category:", category))
        next
      }
      
      BP_gl <- gl_enr[[category]] %>%
        dplyr::rename(adj_pval = Adjusted.P.value) %>%
        dplyr::select(Term, Genes, adj_pval) %>%
        dplyr::mutate(
          Genes = gsub("\\;", ",", Genes),
          Category = category,
          ID = Term,
          regulation = reg
        )
      
      BP_gl <- BP_gl %>%
        dplyr::mutate(
          dbs = category,
          condition = name
        )
      
      regular_terms <- BP_gl %>%
        dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
        dplyr::select(Term, Genes, adj_pval, dbs, condition, regulation) %>%
        tidyr::separate_rows(Genes, sep = ",") %>%
        dplyr::mutate(Genes = trimws(Genes))
      
      regular_terms <- regular_terms %>%
        dplyr::group_by(Term) %>%
        dplyr::filter(dplyr::n_distinct(Genes) >= min_genes_threshold) %>%
        dplyr::ungroup()
      
      top_terms <- BP_gl %>%
        dplyr::rowwise() %>%
        dplyr::mutate(Gene_Count = stringr::str_count(Genes, ",") + 1) %>%
        dplyr::filter(Gene_Count >= min_genes_threshold) %>%
        dplyr::ungroup() %>%
        dplyr::select(-Gene_Count) %>%
        dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
        dplyr::slice_min(adj_pval, n = n) 
      
      filtered_regular <- regular_terms %>%
        dplyr::filter(Term %in% unique(top_terms$Term))
      
      results[[reg]] <- regular_terms
      top_term_list[[reg]]<-top_terms
    }

    regular_results <- dplyr::bind_rows(results[["up-regulated"]], results[["down-regulated"]])
    top_term_results <- dplyr::bind_rows(top_term_list[["up-regulated"]], top_term_list[["down-regulated"]])
    returnlist<-list(regular_results, top_term_results)
    return( returnlist)
  } else {
    gl_enr <- tryCatch(
      enrichr(gl$ID, databases = category, background = background),
      error = function(e) {
        warning(paste("Warning: Enrichment analysis failed for category:", category, "-", e$message))
        return(NULL)
      }
    )
    
    if (is.null(gl_enr) || is.null(gl_enr[[category]]) || nrow(gl_enr[[category]]) == 0) {
      warning(paste("Warning: No enrichment terms found for category:", category))
      return(NULL)
    }
    
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
    
    regular_terms <- BP_gl %>%
      dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
      dplyr::select(Term, Genes, adj_pval, dbs, condition) %>%
      tidyr::separate_rows(Genes, sep = ",") %>%
      dplyr::mutate(Genes = trimws(Genes))
    
    regular_terms <- regular_terms %>%
      dplyr::group_by(Term) %>%
      dplyr::filter(dplyr::n_distinct(Genes) >= min_genes_threshold) %>%
      dplyr::ungroup()
    
    top_terms <- BP_gl %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Gene_Count =  stringr::str_count(Genes, ",") + 1) %>%
      dplyr::filter(Gene_Count >= min_genes_threshold) %>%
      dplyr::ungroup() %>%
      dplyr::select(-Gene_Count) %>%
      dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
      dplyr::slice_min(adj_pval, n = n) 
    
    filtered_regular <- regular_terms %>%
      dplyr::filter(Term %in% unique(top_terms$Term))
    returnlist<-list(regular_terms, top_terms)
    return( returnlist)
  }
}


