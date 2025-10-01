#' Combine P-Values Using Various Methods
#'
#' Pools p-values using one of several methods: Fisher's Method, Stouffer's Method,
#' Weighted Z-test, or Cauchy Combination Test. Caps pooled p-values at a minimum threshold.
#'
#' @param pvallist A numeric vector of p-values.
#' @param method The method to use for combining p-values. Options are `"fisher"`, `"stouffer"`, `"weighted_z"`, and `"cauchy"`.
#' @param weights A numeric vector of weights for the Weighted Z-test. Defaults to equal weights.
#' @param min_pval The minimum threshold for pooled p-values. Defaults to `1e-9`.
#' @return A numeric value representing the combined p-value.
#' @importFrom stats pchisq pnorm qnorm
#' @importFrom VGAM pcauchy
#' @export
poolPValues <- function(pvallist, method = "fisher", weights = NULL, min_pval = 1e-9) {
  # Remove zero p-values (not combinable)
  pvallist <- pvallist[pvallist > 0]  # Avoid log(0) issues
  
  if (length(pvallist) == 0) {
    stop("Error: No valid p-values to combine.")
  }
  
  method <- tolower(method)
  
  # Fisher's Method
  if (method == "fisher") {
    chisq_stat <- -2 * sum(log(pvallist))  # Sum of -2 * log(p-values)
    df <- 2 * length(pvallist)  # Degrees of freedom
    combined_p <- stats::pchisq(chisq_stat, df = df, lower.tail = FALSE)  # Chi-squared test
  } else if (method == "stouffer") {
    # Stouffer's Method
    z_scores <- stats::qnorm(1 - pvallist)
    combined_z <- sum(z_scores) / sqrt(length(z_scores))
    combined_p <- 1 - stats::pnorm(combined_z)
  } else if (method == "weighted_z") {
    # Weighted Z-test
    if (is.null(weights)) {
      weights <- rep(1, length(pvallist))  # Equal weights if not provided
    }
    if (length(weights) != length(pvallist)) {
      stop("Error: weights must be the same length as pvallist.")
    }
    z_scores <- stats::qnorm(1 - pvallist)
    combined_z <- sum(weights * z_scores) / sqrt(sum(weights^2))
    combined_p <- 1 - stats::pnorm(combined_z)
  } else if (method == "cauchy") {
    # Cauchy Combination Test
    t_values <- tan((0.5 - pvallist) * pi)
    combined_cauchy <- mean(t_values)
    combined_p <- VGAM::pcauchy(combined_cauchy, lower.tail = FALSE)
  } else {
    stop("Error: Unsupported method. Choose from 'fisher', 'stouffer', 'weighted_z', or 'cauchy'.")
  }
  
  # Cap pooled p-values at the minimum threshold
  if (combined_p < min_pval) {
    warning(sprintf("Pooled p-value capped to minimum threshold: %.1e", min_pval))
    combined_p <- min_pval
  }
  
  return(combined_p)
}


#' Annotate Clusters with Summary Terms and Pooled P-Values
#'
#' Adds cluster annotations based on a provided term annotation vector and pools p-values
#' for each cluster using the specified method.
#'
#' @param input A data frame containing `condition`, `Cluster`, `Term`, and `adj_pval` columns.
#' @param term_annotation_vector A named vector where names are cluster IDs and values are annotations.
#' @param method The method to use for pooling p-values. Options are `"fisher"`, `"stouffer"`, `"weighted_z"`, and `"cauchy"`.
#' @param weights A numeric vector of weights for the Weighted Z-test. Defaults to `NULL`.
#' @param min_pval The minimum threshold for pooled p-values. Defaults to `1e-9`.
#' @return A data frame with cluster annotations and pooled p-values.
#' @import dplyr
#' @export
annotateClusters <- function(input, term_annotation_vector, method = "fisher", weights = NULL, min_pval = 1e-10) {
 if ("regulation" %in% colnames(input)) {
    annotated_df <- input %>% dplyr::mutate(Cluster_Annotation = term_annotation_vector[as.character(Cluster)]) %>% 
      dplyr::distinct(condition, Cluster, regulation, 
                      Term, .keep_all = TRUE) %>% dplyr::group_by(condition, 
                                                                  Cluster, regulation) %>% dplyr::mutate(pval_pooled = poolPValues(adj_pval, 
                                                                                                                                   method = method, weights = weights, min_pval = min_pval)) %>% 
      dplyr::distinct(condition, Cluster, Term, regulation, 
                      .keep_all = TRUE) %>% ungroup() %>% dplyr::group_by(Cluster, 
                                                                          condition, regulation) %>% dplyr::mutate(unique_terms_per_cluster = list(unique(Term)), 
                                                                                                                         terms_per_cluster = length(unique(Term)), unique_genes_per_cluster = list(unique(unlist(genelist_per_term))), 
                                                                                                                         genes_per_cluster = length(unique(unlist(genelist_per_term))))
    annotated_df <- annotated_df %>%ungroup %>% dplyr::distinct(condition, 
                                                     Cluster, regulation, .keep_all = T) %>% select(-Term, 
                                                                                                    -adj_pval, -dbs, -num_genes_per_term, -genelist_per_term)
  }
  else {
    annotated_df <- input %>% dplyr::mutate(Cluster_Annotation = term_annotation_vector[as.character(Cluster)]) %>% 
      dplyr::distinct(condition, Cluster, 
                      Term, .keep_all = TRUE) %>% dplyr::group_by(condition, 
                                                                  Cluster) %>% dplyr::mutate(pval_pooled = poolPValues(adj_pval, 
                                                                                                                                   method = method, weights = weights, min_pval = min_pval)) %>% 
      dplyr::distinct(condition, Cluster, Term, 
                      .keep_all = TRUE) %>% ungroup() %>% dplyr::group_by(Cluster, 
                                                                          condition) %>% dplyr::mutate(unique_terms_per_cluster = list(unique(Term)), 
                                                                                                                   terms_per_cluster = length(unique(Term)), unique_genes_per_cluster = list(unique(unlist(genelist_per_term))), 
                                                                                                                   genes_per_cluster = length(unique(unlist(genelist_per_term))))
    annotated_df <- annotated_df %>% ungroup %>% dplyr::distinct(condition, 
                                                     Cluster, .keep_all = T) %>% select(-Term, 
                                                                                                    -adj_pval, -dbs, -num_genes_per_term, -genelist_per_term)
  }
  return(annotated_df)
}



#' Plot Enrichment Heatmap
#'
#' Generates a heatmap of enrichment results, with options to split by regulation, customize rotation, colors, 
#' and include additional annotations like terms per cluster.
#'
#' @param input A data frame containing enrichment results with columns `Cluster_Annotation`, `condition`, `pval_pooled`, and optionally `terms_per_cluster`.
#' @param rot Numeric. Rotation angle for column names. Defaults to `0`.
#' @param column_names_centered Logical. If `TRUE`, centers column names. Defaults to `TRUE`.
#' @param cluster_rows Logical. If `TRUE`, rows are clustered. Defaults to `FALSE`.
#' @param cluster_columns Logical. If `TRUE`, columns are clustered. Defaults to `FALSE`.
#' @param legend_name Character. A character for the legend title. Default to "Mean Signif".
#' @param padding simpleUnit. Unt for adjusting margins. 
#' @param plot_colors List. A named list of color palettes for the heatmap. Should contain `up`, `down`, and `default`.
#' @param annotation_bar Logical. If `TRUE`, includes a bar plot annotation for term terms_per_cluster. Defaults to `TRUE`.
#' @return A combined ComplexHeatmap object.
#' @import ComplexHeatmap circlize dplyr tidyr
#' @export
plotHeatmap <- function(
    input,
    rot = 90,
    column_names_centered = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    padding = unit(c(2, 20, 2, 2), "mm"),
    legend_name = "Mean_ Signif",
    plot_colors = list(
      up = c("#ececec", "#e17ecd", "#7900d5"),
      down = c("#ececec", "#41B7C4", "#2A5783"),
      default = c("#ececec", "#41B7C4", "#2A5783")
    ),
    annotation_bar = TRUE
) {
  if("regulation" %in% colnames(input)){split_by_reg = T}
  else{split_by_reg = F}
  #add colum with terms_per_cluster for all conditons
  if(split_by_reg == T){
    terms_per_cluster_all <- input %>%
    dplyr::select(Cluster_Annotation, condition, regulation, unique_terms_per_cluster) %>%
    dplyr::group_by(Cluster_Annotation) %>%
    dplyr::mutate(terms_all = length(unique(unlist(unique_terms_per_cluster)))) %>%
    dplyr::select(-unique_terms_per_cluster)
  # merge it back
  input <- input %>%
    dplyr::left_join(terms_per_cluster_all, by = c("Cluster_Annotation", "condition", "regulation"))
  }
  else{
    terms_per_cluster_all <- input %>%
      dplyr::select(Cluster_Annotation, condition, unique_terms_per_cluster) %>%
      dplyr::group_by(Cluster_Annotation) %>%
      dplyr::mutate(terms_all = length(unique(unlist(unique_terms_per_cluster)))) %>%
      dplyr::select(-unique_terms_per_cluster)
    # merge it back
    input <- input %>%
      dplyr::left_join(terms_per_cluster_all, by = c("Cluster_Annotation", "condition"))    
  }
  prepare_heatmap <- function(data, color_palette, legend_title = legend_name) {
    if (anyNA(data$Cluster_Annotation)) {
      warning("Column 'Cluster_Annotation' contains NA values. Please make sure all clusters are annotated.")
    }
    data_wide <- data %>% tidyr::pivot_wider(names_from = condition, 
                                             values_from = pval_pooled, id_cols = c(Cluster_Annotation, 
                                                                                    terms_all)) %>% dplyr::arrange(desc(terms_all))
    mat <- base::as.matrix(data_wide %>% dplyr::select(-Cluster_Annotation, 
                                                       -terms_all))
    rownames(mat) <- data_wide$Cluster_Annotation
    mat[is.na(mat)] <- 1
    mat[mat == 0] <- 1e-06
    mat <- -log10(mat)
    max_val <- max(mat[is.finite(mat)], na.rm = TRUE)
    mat[is.infinite(mat)] <- max_val + max_val/2
    col_fun <- circlize::colorRamp2(c(0, max_val/2, max_val), 
                                    color_palette)
    row_annotation <- if (annotation_bar && "terms_all" %in% 
                          colnames(data_wide)) {
      ComplexHeatmap::rowAnnotation(NumTerms = ComplexHeatmap::anno_barplot(data_wide$terms_all, 
                                                                            width = grid::unit(4, "cm")))
    }
    else {
      NULL
    }
    col_annotation <- if (annotation_bar && "terms_all" %in% 
                          colnames(data_wide)) {
      ComplexHeatmap::HeatmapAnnotation(NumTerms = ComplexHeatmap::anno_barplot(data_wide$terms_all, 
                                                                                width = grid::unit(4, "cm"), axis_param = list(at = seq(0, 
                                                                                                                                        max(data_wide$terms_all, na.rm = TRUE), 
                                                                                                                                        length.out = 5))), show_annotation_name = split_by_reg && 
                                          legend_title == "Mean Signif down reg. genes")
    }
    else {
      NULL
    }
    if (split_by_reg) {
      heatmap <- ComplexHeatmap::Heatmap(t(mat), name = legend_title, 
                                         col = col_fun, cluster_rows = cluster_rows, 
                                         cluster_columns = cluster_columns, border = FALSE, 
                                         rect_gp = grid::gpar(col = "white", lwd = 1), 
                                         column_names_rot = rot, column_names_centered = column_names_centered, 
                                         column_names_max_height = grid::unit(15, "cm"), 
                                         row_names_max_width = grid::unit(15, "cm"), 
                                         row_names_side = "left")
      if (!is.null(col_annotation)) {
        heatmap <- ComplexHeatmap::Heatmap(t(mat), name = legend_title, 
                                           col = col_fun, cluster_rows = cluster_rows, 
                                           cluster_columns = cluster_columns, border = FALSE, 
                                           rect_gp = grid::gpar(col = "white", lwd = 1), 
                                           column_names_rot = rot, column_names_centered = column_names_centered, 
                                           column_names_max_height = grid::unit(15, "cm"), 
                                           row_names_max_width = grid::unit(15, "cm"), 
                                           row_names_side = "left", top_annotation = col_annotation)
      }
    }
    else {
      heatmap <- ComplexHeatmap::Heatmap(mat, name = "Mean Signif", 
                                         col = col_fun, cluster_rows = cluster_rows, 
                                         cluster_columns = cluster_columns, border = FALSE, 
                                         rect_gp = grid::gpar(col = "white", lwd = 1), 
                                         column_names_rot = rot, column_names_centered = column_names_centered, 
                                         column_names_max_height = grid::unit(15, "cm"), 
                                         row_names_max_width = grid::unit(15, "cm"), 
                                         row_names_side = "left")
      if (!is.null(row_annotation)) {
        heatmap <- heatmap + row_annotation
      }
    }
    return(heatmap)
  }
  if (split_by_reg) {
    data_down <- input %>% dplyr::filter(stringr::str_detect(regulation, 
                                                             "down"))
    data_up <- input %>% dplyr::filter(stringr::str_detect(regulation, 
                                                           "up"))
    all_conditions <- union(unique(data_down$condition), 
                            unique(data_up$condition))
    missing_down <- setdiff(all_conditions, data_down$condition)
    missing_up <- setdiff(all_conditions, data_up$condition)
    if (length(missing_down) > 0) {
      data_down <- dplyr::bind_rows(data_down, tibble::tibble(condition = missing_down, 
                                                              Cluster = data_down$Cluster[1], terms_all = data_down$terms_all[1], 
                                                              Cluster_Annotation = data_down$Cluster_Annotation[1], 
                                                              pval_pooled = NA))
    }
    if (length(missing_up) > 0) {
      data_up <- dplyr::bind_rows(data_up, tibble::tibble(condition = missing_up, 
                                                          Cluster = data_up$Cluster[1], terms_all = data_up$terms_all[1], 
                                                          Cluster_Annotation = data_up$Cluster_Annotation[1], 
                                                          pval_pooled = NA))
    }
    data_down <- data_down %>% dplyr::mutate(condition = factor(condition, 
                                                                levels = all_conditions)) %>% dplyr::arrange(condition)
    data_up <- data_up %>% dplyr::mutate(condition = factor(condition, 
                                                            levels = all_conditions)) %>% dplyr::arrange(condition)
    heatmap_down <- prepare_heatmap(data_down, plot_colors$down, 
                                    legend_title = "Mean Signif down reg. genes")
    heatmap_up <- prepare_heatmap(data_up, plot_colors$up, 
                                  legend_title = "Mean Signif up reg. genes")
    grid::grid.newpage()
    ComplexHeatmap::draw(heatmap_up + heatmap_down, heatmap_legend_side = "right")
    return(invisible(NULL))
  }
  heatmap <- prepare_heatmap(input, plot_colors$default)
  grid::grid.newpage()
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "right")
  return(invisible(NULL))
}


#' Plot Enrichment Bubble plot
#'
#' Generates a Bubble plot of enrichment results, with options to split by regulation
#'
#' @param input A data frame containing enrichment results with columns `Cluster_Annotation`, `condition`, `pval_pooled`, and optionally `terms_per_cluster`.
#' @param plot_colors List. A named list of color palettes for the heatmap. Should contain `up`, `down`, and `default`.
#' @return A ggplot object
#' @import dplyr tidyr ggplot2
#' @export


plotBubbleplot <- function(
    input,
    plot_colors = list(
      up = c("#ececec", "#e17ecd", "#7900d5"),
      down = c("#ececec", "#41B7C4", "#2A5783"),
      default = c("#ececec", "#41B7C4", "#2A5783")
    )
) {
  if ("regulation" %in% colnames(input)) {split_by_reg = T}
  else{split_by_reg = F}
  
  # Prepare data
  input <- input %>%
    dplyr::arrange(pval_pooled) %>%
    dplyr::mutate(Cluster_Annotation = factor(Cluster_Annotation, levels = unique(Cluster_Annotation)))

  # Handle split by regulation
  if (split_by_reg == TRUE) {
    input <- input %>%
      dplyr::mutate(color_condition = ifelse(regulation == "down-regulated", -genes_per_cluster, genes_per_cluster))
    
    Bplot <- ggplot2::ggplot(input, ggplot2::aes(
      x = Cluster_Annotation, y = condition, 
      size = -log10(pval_pooled), color = color_condition
    )) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient2(
        low = plot_colors$down[3], mid = plot_colors$default[1], high = plot_colors$up[3],
        name = "genes per cluster\n(neg. sign = downreg. genes)"
      ) +
      ggplot2::scale_size(name = "-log10(p.adj.)", range = c(3, 10)) + 
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = NULL, title = "") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        legend.position = "left"
      ) +
      ggplot2::facet_wrap(~ regulation)
  }
  else{  
    # Base plot
    Bplot <- ggplot2::ggplot(input, ggplot2::aes(
      x = condition, y = Cluster_Annotation, 
      size = -log10(pval_pooled), color = genes_per_cluster
    )) +
      ggplot2::geom_point() +
      ggplot2::scale_colour_gradient2(low = plot_colors$default[1], mid = plot_colors$default[2], high = plot_colors$default[3]) +
      ggplot2::scale_size(name = "-log10(p.adj.)", range = c(3, 10)) + 
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  
  return(Bplot)
}




