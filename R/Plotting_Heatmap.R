plotting_heatmap <- function(datalist, taxa_filter = NULL, envir_filter = NULL, cluster_cols = F, row_label = NULL,
                         scale_method = make_proportion, col_label = NULL, title = NULL, 
                         col_group = NULL, add_sort = NULL) {
  
  taxa_filter <- enquo(taxa_filter) 
  envir_filter <- enquo(envir_filter)
  col_label <- enquo(col_label)
  row_label <- enquo(row_label)
  col_group <- enquo(col_group)
  add_sort <- enquo(add_sort)
  
  datalist_filtered <- datalist %>%
    with(., if (!rlang::quo_is_null(taxa_filter)) filter_taxa_datalist(., !!taxa_filter) else .) %>%
    with(., if (!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .) %>%
    with(., suppressWarnings(arrange_datalist(., !!col_group, !!add_sort, !!col_label)))
  
  if (nrow(datalist_filtered$Count_Data) == 0) {stop("No taxa with such name.")}
  
  if (!rlang::quo_is_null(col_group)) {
    if (nrow(unique(select(datalist_filtered$Meta_Data, !!col_group))) > 1) {
      gaps_col <- cumsum(table(select(datalist_filtered$Meta_Data, !!col_group, !!add_sort)))
    } else {
      gaps_col <- NULL
    }
  } else {
    gaps_col <- NULL
  }
  
  Count_Matrix <- datalist_filtered %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    as.matrix() %>%
    apply(., 1, function(x) scale_method(x)) %>% 
    t()
  
  rownames(Count_Matrix) <- datalist_filtered$Count_Data %>%
    select(1) %>%
    deframe()
  
  row_label <- if(!rlang::quo_is_null(row_label)) {
    datalist_filtered$Count_Data %>%
      select(!!row_label) %>%
      deframe()
  } else NULL
  
  col_label <- if (!rlang::quo_is_null(col_label)) {
    
      datalist_filtered$Meta_Data %>%
      #  arrange(!!col_group, !!col_label) %>%
        mutate(!!col_label := if (is.numeric(!!col_label)) round(!!col_label, digits = 0) else !!col_label) %>%
        select(!!col_label) %>%
        deframe() %>%
        as.character()
    
  } else NULL
  
  if (nrow(Count_Matrix) > 2) {
    clust_num <- factoextra::fviz_nbclust(Count_Matrix, cluster::pam, method = "silhouette", k.max = max(min(nrow(Count_Matrix)-1, 8),2))$data$y %>%
      which.max()
    
    clust_num <- 4
    
    cluster_obj <- Count_Matrix %>%
      vegan::vegdist(method = "euclidian") %>%
      hclust(., method = "ward.D2") %>%
      with(., if ("Latitude" %in% names(datalist_filtered$Meta_Data)) {
        reorder(., wts = abs(datalist_filtered$Meta_Data$Latitude), agglo.FUN = "uwmean")
      } else .)
    
  } else {
    
    cluster_obj <- FALSE
    
  }
  
  row_annotation <- data.frame(Cluster = paste("Cluster", if (is.logical(cluster_obj)) rep(1, nrow(Count_Matrix)) else cutree(cluster_obj, k = clust_num)), row.names = rownames(Count_Matrix))
  row_annotation_color <- list(Cluster = structure(ggsci::pal_rickandmorty(palette = c("schwifty"))(ifelse(is.logical(cluster_obj), 1, clust_num)), names = unique(row_annotation$Cluster)))
  
  if (is.null(title)) {
    title <- paste(str_split_fixed(as_label(taxa_filter), pattern = "==", 2)[1] %>% str_replace(., " ", ""), ": ",
                   str_split_fixed(as_label(taxa_filter), pattern = "==", 2)[2] %>% str_replace_all(., "\\\"", "") %>% str_replace(., " ", ""), sep = "")
  }
  
  pheatmap(Count_Matrix, cluster_cols = cluster_cols, cluster_rows = cluster_obj, cutree_rows = clust_num, #clustering_method = "ward.D2", 
           main = title, border_color = "grey60",
           legend_breaks = c(min(Count_Matrix), max(Count_Matrix)),legend_labels = c("Min  ", "Max  "),
           angle_col = "90", annotation_legend = F, labels_col = col_label,
           labels_row = row_label, show_rownames = ifelse(!is.null(row_label),T,F),
           annotation_row = row_annotation, fontsize_row = 6, fontsize_col = 7,
           annotation_colors = row_annotation_color, fontsize = 7,
           gaps_col = gaps_col)
  
}