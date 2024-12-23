#' Evaluate Thresholds for Network Clustering
#'
#' Evaluates graph properties (connected components, modularity, cluster count) for a range of thresholds.
#'
#' @param input A data frame containing terms, genes, and additional metadata.
#' @return A Plotly object visualizing the evaluation of thresholds.
#' @import dplyr tidyr igraph plotly proxy
#' @export
evaluateThreshold <- function(input, seed = 42) {
  # Ensure necessary columns exist
  set.seed(seed)
  required_columns <- c("Term", "Genes")
  missing_columns <- base::setdiff(required_columns, colnames(input))
  if (length(missing_columns) > 0) {
    stop("Input data is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  
  # Preprocessing
  input <- input %>% dplyr::select(-adj_pval, everything())
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

  # Define thresholds and initialize results
  thresholds <- base::seq(0.05, 1, by = 0.01)
  results <- base::data.frame(
    Threshold = thresholds,
    Connected_Components = base::numeric(length(thresholds)),
    Modularity = base::numeric(length(thresholds)),
    Cluster_Count = base::numeric(length(thresholds))
  )
  
  # Evaluate metrics for each threshold
  for (i in base::seq_along(thresholds)) {
    threshold <- thresholds[i]
    
    # Create graph
    graph <- igraph::graph_from_adjacency_matrix(base::as.matrix(distance_matrix), mode = "undirected", weighted = TRUE)
    graph <- igraph::delete_edges(graph, igraph::E(graph)[igraph::E(graph)$weight < threshold])
    
    # Check for valid graph
    if (igraph::vcount(graph) > 0 && igraph::ecount(graph) > 0) {
      components_result <- igraph::components(graph)
      clusters <- igraph::cluster_walktrap(graph)
      results$Connected_Components[i] <- components_result$no
      results$Modularity[i] <- igraph::modularity(clusters)
      results$Cluster_Count[i] <- length(clusters)
    } else {
      results$Connected_Components[i] <- 0
      results$Modularity[i] <- NA
      results$Cluster_Count[i] <- 0
    }
  }
  
  # Plot results
  results_plot <- results %>%
    tidyr::pivot_longer(cols = -Threshold, names_to = "Metric", values_to = "Value") %>%
    plotly::plot_ly(
      x = ~Threshold, y = ~Value, color = ~Metric,
      colors = c("#00B2AE", "#6247A2", "#E63946")
    ) %>%
    plotly::add_lines() %>%
    plotly::layout(
      title = "Threshold Selection Plot",
      xaxis = list(title = "Threshold"),
      yaxis = list(title = "Metric Value"),
      legend = list(title = list(text = "Metric")),
      template = "plotly_white"
    )
  
  return(results_plot)
}


