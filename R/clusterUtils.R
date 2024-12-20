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
  annotated_df <- input %>%
    dplyr::mutate(Cluster_Annotation = term_annotation_vector[as.character(Cluster)]) %>%
    dplyr::distinct(condition, Cluster, Term, .keep_all = TRUE) %>%
    dplyr::group_by(condition, Cluster) %>%
    dplyr::mutate(pval_pooled = poolPValues(adj_pval, method = method, weights = weights, min_pval = min_pval)) %>%
    dplyr::distinct(condition, Cluster_Annotation, pval_pooled, .keep_all = TRUE)
  
  return(annotated_df)
}


#' Plot Enrichment Heatmap
#'
#' Generates a heatmap of enrichment results, with options to split by regulation, customize rotation, colors, 
#' and include additional annotations like terms per cluster.
#'
#' @param input A data frame containing enrichment results with columns `Cluster_Annotation`, `condition`, `pval_pooled`, and optionally `terms_per_cluster`.
#' @param split_by_reg Logical. If `TRUE`, splits heatmaps into up-regulated and down-regulated genes.
#' @param rot Numeric. Rotation angle for column names. Defaults to `0`.
#' @param column_names_centered Logical. If `TRUE`, centers column names. Defaults to `TRUE`.
#' @param cluster_rows Logical. If `TRUE`, rows are clustered. Defaults to `FALSE`.
#' @param cluster_columns Logical. If `TRUE`, columns are clustered. Defaults to `FALSE`.
#' @param legend_name Character. A character for the legend title. Default to "Mean Signif".
#' @param colors List. A named list of color palettes for the heatmap. Should contain `up`, `down`, and `default`.
#' @param annotation_bar Logical. If `TRUE`, includes a bar plot annotation for term terms_per_cluster. Defaults to `TRUE`.
#' @return A combined ComplexHeatmap object.
#' @import ComplexHeatmap circlize dplyr tidyr
#' @export
plotHeatmap <- function(
    input,
    split_by_reg = FALSE,
    rot = 0,
    column_names_centered = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    legend_name = "Mean_ Signif",
    plot_colors = list(
      up = c("#ececec", "#e17ecd", "#7900d5"),
      down = c("#ececec", "#41B7C4", "#2A5783"),
      default = c("#ececec", "#41B7C4", "#2A5783")
    ),
    annotation_bar = TRUE
) {
  # Helper function to prepare data and generate a heatmap
  prepare_heatmap <- function(data, color_palette, legend_title = legend_name) {
    # Pivot data to wide format
    data_wide <- data %>%
      tidyr::pivot_wider(names_from = condition, values_from = pval_pooled, id_cols = c(Cluster_Annotation, terms_per_cluster)) %>%
      dplyr::arrange(desc(terms_per_cluster))  # Sort rows by terms_per_cluster
    
    # Create matrix for heatmap
    mat <- base::as.matrix(data_wide %>% dplyr::select(-Cluster_Annotation, -terms_per_cluster))
    rownames(mat) <- data_wide$Cluster_Annotation
    mat[is.na(mat)] <- 1  # Fill missing values with p-value = 1
    mat[mat == 0] <- 1e-6  # Avoid log(0)
    mat <- -log10(mat)
    max_val <- max(mat[is.finite(mat)], na.rm = TRUE)
    mat[is.infinite(mat)] <- max_val + max_val / 2  # Cap infinite values
    
    # Create color function
    col_fun <- circlize::colorRamp2(c(0, max_val / 2, max_val), color_palette)
    
    # Add optional frequency annotation
    row_annotation <- if (annotation_bar && "terms_per_cluster" %in% colnames(data_wide)) {
      ComplexHeatmap::rowAnnotation(NumTerms = ComplexHeatmap::anno_barplot(data_wide$terms_per_cluster, width = grid::unit(4, "cm")))
      
    } else {
      NULL
    }
    col_annotation <- if (annotation_bar && "terms_per_cluster" %in% colnames(data_wide)) {

      ComplexHeatmap::HeatmapAnnotation(
        NumTerms = ComplexHeatmap::anno_barplot(
          data_wide$terms_per_cluster,
          width = grid::unit(4, "cm"),
          axis_param = list(at = seq(0, max(data_wide$terms_per_cluster, na.rm = TRUE), length.out = 5))
        ),
        show_annotation_name = split_by_reg && legend_title == "Mean Signif down reg. genes"
      )
    } else {
      NULL
    }
    # Generate heatmap
    heatmap <- ComplexHeatmap::Heatmap(
      mat, name = "Mean Signif", col = col_fun, cluster_rows = cluster_rows, cluster_columns = cluster_columns,
      border = FALSE, rect_gp = grid::gpar(col = "white", lwd = 1), column_names_rot = rot,
      column_names_centered = column_names_centered, column_names_max_height = grid::unit(15, "cm"),
      row_names_max_width = grid::unit(15, "cm"), row_names_side = "left"
    )
    # Combine heatmap with row annotation if applicable
    if (!is.null(row_annotation)) {
      heatmap <- heatmap + row_annotation
    }
    
    # Heatap must be turned to combine up- and downregulation
    if (!is.null(split_by_reg)) {
      heatmap <- ComplexHeatmap::Heatmap(
        t(mat), name = legend_title, col = col_fun, cluster_rows = cluster_rows, cluster_columns = cluster_columns,
        border = FALSE, rect_gp = grid::gpar(col = "white", lwd = 1), column_names_rot = rot,
        column_names_centered = column_names_centered, column_names_max_height = grid::unit(15, "cm"),
        row_names_max_width = grid::unit(15, "cm"), row_names_side = "left"
      )
      if (!is.null(col_annotation)) {
        heatmap <- ComplexHeatmap::Heatmap(
          t(mat), name = legend_title, col = col_fun, cluster_rows = cluster_rows, cluster_columns = cluster_columns,
          border = FALSE, rect_gp = grid::gpar(col = "white", lwd = 1), column_names_rot = rot,
          column_names_centered = column_names_centered, column_names_max_height = grid::unit(15, "cm"),
          row_names_max_width = grid::unit(15, "cm"), row_names_side = "left",  top_annotation = col_annotation)
      }
      
    }
    
    
    
    return(heatmap)
  }
  
  if (split_by_reg) {
    # Split input into up-regulated and down-regulated data
    data_down <- input %>%
      dplyr::filter(stringr::str_detect(regulation, "down")) 
    
    data_up <- input %>%
      dplyr::filter(stringr::str_detect(regulation, "up")) 
    
    # Ensure both datasets have matching conditions
    all_conditions <- union(unique(data_down$condition), unique(data_up$condition))
    missing_down <- setdiff(all_conditions, data_down$condition)
    missing_up <- setdiff(all_conditions, data_up$condition)
    
    if (length(missing_down) > 0) {
      data_down <- dplyr::bind_rows(data_down, tibble::tibble(condition = missing_down, Cluster = NA, Cluster_Annotation = NA, pval_pooled = NA))
    }
    if (length(missing_up) > 0) {
      data_up <- dplyr::bind_rows(data_up, tibble::tibble(condition = missing_up, Cluster = NA, Cluster_Annotation = NA, pval_pooled = NA))
    }
    
    # Reorder both datasets to match the order of all_conditions
    data_down <- data_down %>%
      dplyr::mutate(condition = factor(condition, levels = all_conditions)) %>%
      dplyr::arrange(condition)
    
    data_up <- data_up %>%
      dplyr::mutate(condition = factor(condition, levels = all_conditions)) %>%
      dplyr::arrange(condition)
    
    # Generate heatmaps for up and down regulation
    heatmap_down <- prepare_heatmap(data_down, plot_colors$down, legend_title = "Mean Signif down reg. genes")
    heatmap_up <- prepare_heatmap(data_up, plot_colors$up,  legend_title = "Mean Signif up reg. genes")
    
    
    # Combine heatmaps and force plot
    grid::grid.newpage()
    ComplexHeatmap::draw(heatmap_up + heatmap_down, heatmap_legend_side = "right")
    return(invisible(NULL))
  }
  
  # Single heatmap for non-split data
  heatmap <- prepare_heatmap(input, plot_colors$default)
  grid::grid.newpage()
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "right")  # Explicitly draw heatmap
  return(invisible(NULL))
}


