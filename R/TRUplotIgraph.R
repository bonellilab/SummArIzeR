#' Generate an Interactive Network Graph with Repelled Labels
#'
#' Creates a Plotly-based interactive graph visualization of terms and their relationships, minimizing label overlap.
#'
#' @param input A data frame containing terms, genes, and additional metadata.
#' @param ts A numeric threshold for edge inclusion based on distance weights.
#' @return A Plotly object representing the network graph.
#' @import dplyr tidyr igraph plotly proxy factoextra
#' @export
TRUplotIgraph <- function(input, ts, seed = 42) {
  # Preprocessing
  set.seed(seed)
  input <- input %>% dplyr::select(-adj_pval)
  gene_freq <- as.data.frame(input %>%
                               dplyr::group_by(Term,  Genes) %>%
                               dplyr::summarise(Frequency = dplyr::n(), .groups = "drop") %>%
                               tidyr::pivot_wider(names_from = c(Term), values_from = Frequency, values_fill = 0))
  
  rownames(gene_freq) <- gene_freq$Genes
  gene_freq <- gene_freq %>% dplyr::select(-Genes)
  # Create binary matrix
  matrix_genes <- as.matrix(gene_freq)
  matrix_genes[matrix_genes >= 1] <- 1

  # Generate similarity and distance matrices
  distance_matrix <- proxy::dist(t(matrix_genes), method = "Jaccard", convert_similarities = F)
  print(fviz_dist(
    distance_matrix,
    order = TRUE,
    show_labels = TRUE,
    lab_size = NULL,
    gradient = list(low = "red", mid = "white", high = "blue")
  ))
  # Create graph
  graph <- igraph::graph_from_adjacency_matrix(as.matrix(distance_matrix), mode = "undirected", weighted = TRUE)
  graph <- igraph::delete_edges(graph, igraph::E(graph)[weight < ts])
  
  clusters <- igraph::cluster_walktrap(graph, merges = T)
  igraph::plot_dendrogram(clusters)
  # Generate color palette
  palette <- rainbow(max(igraph::membership(clusters)))
  igraph::V(graph)$color <- palette[igraph::membership(clusters)]
  
  # Extract graph layout
  layout <- igraph::layout_with_fr(graph)  # Force-directed layout
  
  # Convert layout to data frame for Plotly
  layout_df <- as.data.frame(layout)
  colnames(layout_df) <- c("x", "y")
  layout_df$node <-  igraph::V(graph)$name
  layout_df$color <-  igraph::V(graph)$color
  layout_df$cluster <- igraph::membership(clusters)[match(layout_df$node,  igraph::V(graph)$name)]
  # Calculate number of genes per term
  term_gene_counts <- input %>%
    dplyr::group_by(Term) %>%
    dplyr::summarise(num_genes = dplyr::n_distinct(Genes), .groups = "drop")
  
  layout_df <- layout_df %>%
    dplyr::left_join(term_gene_counts, by = c("node" = "Term"))
 
   # Insert database information
  dbs_info_df <- input %>%
    dplyr::group_by(Term) %>%
    dplyr::reframe(dbs_info = dbs)
  
  
  layout_df <- layout_df %>%
    dplyr::left_join(dbs_info_df, by = c("node" = "Term"))
  
  # Add hover information
  layout_df$hover_info <- paste0(
    "Term: ", layout_df$node, "<br>",
    "Database: ", layout_df$dbs_info, "<br>",
    "Cluster: ", layout_df$cluster, "<br>",
    "Genes: ", layout_df$num_genes
  )
  
  # Create edges for Plotly
  edges <- igraph::as_data_frame(graph, what = "edges")
  edges <- edges %>%
    dplyr::mutate(
      x_start = layout_df$x[match(from, layout_df$node)],
      y_start = layout_df$y[match(from, layout_df$node)],
      x_end = layout_df$x[match(to, layout_df$node)],
      y_end = layout_df$y[match(to, layout_df$node)]
    )
  
  # Plot with Plotly
  plot <- plotly::plot_ly() %>%
    # Add edges
    plotly::add_segments(
      data = edges,
      x = ~x_start, y = ~y_start, xend = ~x_end, yend = ~y_end,
      line = list(color = 'gray', width = 0.5),
      showlegend = FALSE
    ) %>%
    # Add nodes with hover text
    plotly::add_trace(
      data = layout_df,
      x = ~x, y = ~y,
      type = "scatter",
      mode = "markers",
      marker = list(size = 10, color = ~color),
      text = ~hover_info,  # Hover text with term, cluster, and gene count
      hoverinfo = "text",  # Display hover text only
      showlegend = FALSE
    ) %>%
    plotly::layout(
      title = "Interactive Network Graph",
      xaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", zeroline = FALSE, showgrid = FALSE, showticklabels = FALSE)
    )
  
  return(plot)
}
