#' Assign Clusters to Terms Using Graph Clustering
#'
#' Performs graph-based clustering and assigns clusters to terms, including metadata and summary statistics.
#'
#' @param input A data frame containing terms, genes, and additional metadata.
#' @param ts A numeric threshold for edge inclusion based on distance weights.
#' @return A data frame with terms, clusters, and enriched metadata, including a summary of terms per cluster.
#' @import dplyr tidyr igraph
#' @export
returnIgraphCluster<- function(input, ts, seed = 42) {


  # Ensure necessary columns exist
  set.seed(seed)
  required_columns <- c("Term", "Genes", "dbs")
  missing_columns <- base::setdiff(required_columns, colnames(input))
  if (length(missing_columns) > 0) {
    stop("Input data is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # Preprocessing
  gene_freq <- input %>%
    dplyr::group_by(Term, Genes) %>%
    dplyr::summarise(Frequency = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Term, values_from = Frequency, values_fill = 0)
  
  gene_freq <- base::as.data.frame(gene_freq)  # Convert tibble to data frame
  matrix_genes <- base::as.matrix(gene_freq %>% dplyr::select(-Genes))  # Binary matrix
  rownames(matrix_genes) <- gene_freq$Genes  # Assign row names
  matrix_genes[matrix_genes >= 1] <- 1  # Binary transformation
  
  # Generate similarity and distance matrices
  distance_matrix <- proxy::dist(t(matrix_genes), method = "Jaccard", convert_similarities = F)

  # Create graph
  graph <- igraph::graph_from_adjacency_matrix(base::as.matrix(distance_matrix), mode = "undirected", weighted = TRUE)
  graph <- igraph::delete_edges(graph, igraph::E(graph)[igraph::E(graph)$weight < ts])
  
  # Cluster analysis
  clusters <- igraph::cluster_walktrap(graph)
  cluster_membership <- igraph::membership(clusters)
  
  # Add cluster information to the terms
  term_cluster <- data.frame(
    Term = igraph::V(graph)$name,
    Cluster = as.vector(cluster_membership)
  )
  
  # Add metadata:`dbs`
  term_cluster <- term_cluster %>%
    dplyr::mutate(
      dbs = sapply(Term, function(term_name) {
        input %>%
          dplyr::filter(Term == term_name) %>%
          dplyr::pull(dbs) %>%
          base::unique() %>%
          base::paste(collapse = ", ")
      })
    )


  Termlist_all <- input %>%
    dplyr::left_join(term_cluster, by = "Term") %>%  # Merge clusters
    dplyr::select(-dbs.x) %>%  # Remove redundant dbs column
    dplyr::rename(dbs = dbs.y) %>%  # Rename dbs column
    dplyr::group_by(condition, Term) %>%  
    dplyr::mutate(
      num_genes_per_term = length(unique(unlist(Genes))),  # Count unique genes
      genelist_per_term = paste(unique(unlist(na.omit(Genes))), collapse = ", "),  # Concatenate unique genes
    ) %>%
  dplyr::distinct(condition, Cluster,Term,adj_pval, .keep_all = T ) %>% #Keep distinct Terms for every conditions
    dplyr::select(-Genes)  
  
  return(Termlist_all)
}



