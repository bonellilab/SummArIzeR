#' Generate ChatGPT Prompts for Cluster Summarization
#'
#' Creates ChatGPT prompts for summarizing terms in each cluster and provides instructions
#' for assigning cluster numbers to summary terms.
#'
#' @param input A data frame containing `Term` and `Cluster` columns.
#' @return None. Outputs ChatGPT prompts to the console.
#' @import dplyr
#' @export
generateGPTPrompt <- function(input) {
  # Extract unique terms per cluster
  unique_terms_per_cluster <- input %>%
    dplyr::ungroup() %>%
    dplyr::select(Term, Cluster) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Cluster)
  
  # Generate ChatGPT prompts for each cluster
  clusters <- unique_terms_per_cluster %>%
    dplyr::pull(Cluster) %>%
    unique()
  
  for (cluster in clusters) {
    terms <- unique_terms_per_cluster %>%
      dplyr::filter(Cluster == cluster) %>%
      dplyr::pull(Term)
    terms_string <- paste(terms, collapse = ", ")
    cat(sprintf("Cluster %d:\nPlease find a summary term for the following terms: %s\n\n", cluster, terms_string))
  }
  
  # Print final instruction for ChatGPT
  cat("Please summarize the following clusters of terms with a single term for each cluster. \n",
      "Then, assign the cluster numbers to these summary terms in an R vector, called `cluster_summary`, formatted as:\n",
      "'1' = 'Summary Term 1', '2' = 'Summary Term 2', ..., 'n' = 'Summary Term n'.\n")
}
