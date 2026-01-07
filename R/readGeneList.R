# ----------------------------
# Internal GMT enrichment helpers
# ----------------------------

#' @keywords internal
.read_gmt_pathways <- function(gmt_files) {
  gmt_files <- as.character(gmt_files)
  stopifnot(length(gmt_files) >= 1)
  
  gsets_list <- lapply(gmt_files, function(f) {
    if (requireNamespace("fgsea", quietly = TRUE)) {
      fgsea::gmtPathways(f)
    } else {
      lines <- readLines(f, warn = FALSE)
      sets <- lapply(lines, function(x) {
        parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
        genes <- unique(parts[-c(1, 2)])
        genes[nzchar(genes)]
      })
      names(sets) <- vapply(
        lines,
        function(x) strsplit(x, "\t", fixed = TRUE)[[1]][1],
        character(1)
      )
      sets
    }
  })
  
  names(gsets_list) <- sub("\\.gmt$", "", basename(gmt_files), ignore.case = TRUE)
  gsets_list
}

#' @keywords internal
.enrich_gmt_like_enrichr <- function(
    genes,
    gmt_files,
    background = NULL,
    min_overlap = 1,
    min_set_size = 5,
    max_set_size = 5000,
    default_universe_size = 20000,
    count_unmapped_in_input = TRUE,   # Enrichr-like behavior
    verbose = FALSE,
    show_progress = interactive()
) {
  genes <- unique(as.character(genes))
  genes <- genes[nzchar(genes)]
  if (length(genes) == 0) return(data.frame())
  
  # Read pathways (list: db -> list(term -> genes))
  gsets_by_file <- .read_gmt_pathways(gmt_files)
  
  # Flatten into parallel vectors
  db <- rep(names(gsets_by_file), times = vapply(gsets_by_file, length, integer(1)))
  term <- unlist(lapply(gsets_by_file, names), use.names = FALSE)
  sets <- unlist(gsets_by_file, recursive = FALSE, use.names = FALSE)
  
  # Mapping universe for overlaps and set membership
  if (is.null(background)) {
    background_genes <- unique(unlist(sets, use.names = FALSE))
    background_genes <- background_genes[nzchar(background_genes)]
    N <- as.integer(default_universe_size)
  } else {
    background_genes <- unique(as.character(background))
    background_genes <- background_genes[nzchar(background_genes)]
    N <- length(background_genes)
  }
  
  # Mapped genes (used for overlaps)
  genes_bg <- intersect(genes, background_genes)
  
  # Enrichr-like: count input genes even if they don't map into the database
  if (count_unmapped_in_input) {
    n_input <- length(genes)
    n_input <- min(n_input, N)  # hypergeometric requires n_input <= N
  } else {
    n_input <- length(genes_bg)
  }
  
  # Prepare sets restricted to mapping universe
  sets_bg <- lapply(sets, function(s) intersect(unique(as.character(s)), background_genes))
  set_sizes <- vapply(sets_bg, length, integer(1))
  
  keep <- set_sizes >= min_set_size & set_sizes <= max_set_size
  if (!any(keep)) return(data.frame())
  
  db <- db[keep]
  term <- term[keep]
  sets_bg <- sets_bg[keep]
  set_sizes <- set_sizes[keep]
  
  # Safety: hypergeometric requires m <= N and n_input <= N
  # If your GMT implies m > N, bump N so the model is coherent.
  if (is.null(background)) {
    max_m <- max(set_sizes)
    if (max_m > N) {
      if (verbose) message("GMT enrichment: increasing universe N from ", N, " to ", max_m,
                           " because a gene set size exceeds N.")
      N <- as.integer(max_m)
      n_input <- min(n_input, N)
    }
  }
  
  total_tests <- length(sets_bg)
  if (verbose) {
    message(
      "GMT enrichment: ",
      length(unique(db)), " database(s), ",
      total_tests, " gene sets | ",
      "input=", length(genes),
      " mapped=", length(genes_bg),
      " n_input=", n_input,
      " N=", N,
      if (is.null(background)) " (default universe)" else " (provided background)"
    )
  }
  
  pb <- NULL
  if (show_progress) {
    pb <- utils::txtProgressBar(min = 0, max = total_tests, style = 3)
    on.exit(try(close(pb), silent = TRUE), add = TRUE)
  }
  
  out <- vector("list", total_tests)
  
  for (i in seq_along(sets_bg)) {
    if (show_progress) utils::setTxtProgressBar(pb, i)
    
    gs <- sets_bg[[i]]
    m <- set_sizes[[i]]
    
    # overlap only among mappable genes
    a <- sum(genes_bg %in% gs)
    if (a < min_overlap) next
    
    # one-sided enrichment p-value: P[X >= a], X ~ Hypergeom(N, m, n_input)
    p <- stats::phyper(a - 1, m, N - m, n_input, lower.tail = FALSE)
    
    # odds ratio (Enrichr-style) uses the 2x2 table:
    # a overlap
    # b in set only
    # c in input only
    # d neither
    #
    # If count_unmapped_in_input=TRUE, c includes unmapped input genes as well.
    b <- m - a
    c <- n_input - a
    d <- N - (a + b + c)
    if (d < 0) d <- 0  # guard against edge cases
    
    odds_ratio <- (1.0 * a * d) / max(1.0 * b * c, 1)
    combined <- -log(p) * odds_ratio
    
    overlap_genes <- intersect(genes_bg, gs)
    
    out[[i]] <- data.frame(
      Database = db[i],
      Term = term[i],
      P.value = p,
      Adjusted.P.value = NA_real_,
      Odds.Ratio = odds_ratio,
      Combined.Score = combined,
      Genes = paste(overlap_genes, collapse = ";"),
      Overlap = paste0(a, "/", m),
      stringsAsFactors = FALSE
    )
  }
  
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(data.frame())
  
  # Adjust p-values within each database then merge
  out$Adjusted.P.value <- ave(
    out$P.value,
    out$Database,
    FUN = function(p) stats::p.adjust(p, method = "BH")
  )
  
  out <- out[order(out$Database, out$Adjusted.P.value, -out$Combined.Score), , drop = FALSE]
  rownames(out) <- NULL
  
  if (verbose) message("GMT enrichment complete.")
  out
}



#' Read Gene List and Perform Enrichment
#'
#' Reads a gene list, optionally splits by regulation, and performs enrichment analysis.
#'
#' @param listpathSig A data frame containing the gene list with columns `genes` and optionally `log2fold`.
#' @param name The name of the condition (used in results).
#' @param category The enrichment database to use.
#' @param background A background list of genes to use with enrichR. See enrichr documentation for details. 
#' @param split_by_reg Logical. If 'TRUE', splits genes by up- and down-regulation.
#' @param logFC_threshold The threshold for determining up/down-regulation. Defaults to `0`.
#' @param pval_threshold The adjusted p-value threshold for filtering results.
#' @param n The number of top terms to return. Default is 10.
#' @param min_genes_threshold The minimum number of genes required for a term to be kept in the results. Default is 3. 
#' @param enrichment_method Which method to use for enrichment, "enrichr" (requires internet connection), or "gmt" (offline, requires gmt file input. 
#' @param gmt_files Required if enrichment_method = "gmt". Path to one or more local gmt file(s). For gmt data format see https://docs.gsea-msigdb.org/#GSEA/Data_Formats/. 
#' @param gmt_min_overlap Minimum overlap between geneset and dataset
#' @param gmt_min_set_size Minimum size of geneset
#' @param gmt_max_set_size Maximum size of geneset
#' @param verbose Logical. If TRUE, print informative messages about the progress. 
#' @param show_progress Logical. If TRUE, display a progress bar during computation.
#' @return A list with two data frames: `regular` (all terms) and `top` (top terms).
#' @import dplyr tidyr
#' @import stringr
#' @export

readGeneList <- function(listpathSig, name, category, background = NULL, split_by_reg = FALSE, logFC_threshold = 0,
                         pval_threshold = 0.05, n = 10, min_genes_threshold = 3,
                         enrichment_method = c("enrichr", "gmt"),
                         gmt_files = NULL,
                         gmt_min_overlap = 3,
                         gmt_min_set_size = 5,
                         gmt_max_set_size = 5000,
                         verbose = FALSE,
                         show_progress = FALSE) {
  
  enrichment_method <- match.arg(enrichment_method)
  
  if (nrow(listpathSig) == 0) stop("Error: The input gene list (listpathSig) is empty.")
  if (!"genes" %in% colnames(listpathSig)) stop("Error: The input gene list must contain a 'genes' column.")
  
  gl <- listpathSig %>% dplyr::rename("ID" = "genes")
  
  if (verbose) {
    msg_reg <- if (split_by_reg) "split" else "unsplit"
    message("readGeneList: method=", enrichment_method,
            " | condition=", name,
            " | category=", category,
            " | genes=", nrow(gl),
            " | ", msg_reg)
  }
  
  .run_enrichment <- function(ids, category, background) {
    
    if (enrichment_method == "enrichr") {
      tryCatch(
        enrichr(ids, databases = category, background = background),
        error = function(e) {
          warning(paste("Warning: Enrichment analysis failed for category:", category, "-", e$message))
          return(NULL)
        }
      )
    } else {
      if (is.null(gmt_files) || length(gmt_files) == 0) {
        stop("Error: enrichment_method='gmt' requires 'gmt_files' (a single GMT file after mapping).")
      }
      
      df <- .enrich_gmt_like_enrichr(
        genes = ids,
        gmt_files = gmt_files,  # already per-category
        background = background,
        min_overlap = gmt_min_overlap,
        min_set_size = gmt_min_set_size,
        max_set_size = gmt_max_set_size,
        verbose = verbose,
        show_progress = show_progress
      )
      
      if (is.null(df) || nrow(df) == 0) return(NULL)
      
      # Mimic enrichR output: list keyed by category label
      out <- df
      out$Genes <- gsub("\\;", ",", out$Genes)
      
      names_list <- list(out)
      names(names_list) <- category
      names_list
    }
  }
  
  if (split_by_reg) {
    if (!"log2fold" %in% colnames(gl)) stop("Error: 'log2fold' column is required for splitting by regulation.")
    
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
    top_term_list <- list()
    
    for (reg in c("up-regulated", "down-regulated")) {
      gl_subset <- gl %>% dplyr::filter(regulation == reg)
      if (verbose) message("  readGeneList: running ", reg, " (n=", nrow(gl_subset), ")")
      
      gl_enr <- .run_enrichment(gl_subset$ID, category = category, background = background)
      
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
        ) %>%
        dplyr::mutate(dbs = category, condition = name)
      
      regular_terms <- BP_gl %>%
        dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
        dplyr::select(Term, Genes, adj_pval, dbs, condition, regulation) %>%
        tidyr::separate_rows(Genes, sep = ",") %>%
        dplyr::mutate(Genes = trimws(Genes)) %>%
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
      
      results[[reg]] <- regular_terms
      top_term_list[[reg]] <- top_terms
    }
    
    regular_results <- dplyr::bind_rows(results[["up-regulated"]], results[["down-regulated"]])
    top_term_results <- dplyr::bind_rows(top_term_list[["up-regulated"]], top_term_list[["down-regulated"]])
    return(list(regular_results, top_term_results))
    
  } else {
    
    gl_enr <- .run_enrichment(gl$ID, category = category, background = background)
    
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
      ) %>%
      dplyr::mutate(dbs = category, condition = name)
    
    regular_terms <- BP_gl %>%
      dplyr::filter(!is.na(adj_pval) & adj_pval < pval_threshold) %>%
      dplyr::select(Term, Genes, adj_pval, dbs, condition) %>%
      tidyr::separate_rows(Genes, sep = ",") %>%
      dplyr::mutate(Genes = trimws(Genes)) %>%
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
    
    return(list(regular_terms, top_terms))
  }
}
