
weighted_mean_datalist <- function(datalist, param) {
  param <- enquo(param)
  otu_id <- pull(datalist$Count_Data, 1)
  param_val <- pull(datalist$Meta_Data, !!param)
  
  final <- map(otu_id, function(x) {
    tmp <- weighted.mean(x = param_val, w = filter(datalist$Count_Data, pull(datalist$Count_Data, 1) == !!x) %>%
                    select_if(is.numeric) %>%
                    as.numeric(), na.rm = T)
    return(tibble(Cluster = rlang::as_name(x), !!rlang::quo_name(param) := tmp))
  }) %>%
    bind_rows()
  
  return(final)
}