sparCC_wrapper <- function(datalist, envir_filter = NULL, n_boot = 99, frac = FALSE) {
  
  shuffle <- function(count) {
    
    tmp <- count
    
    for (i in 1:ncol(tmp)) {
      tmp[,i] <- tmp[sample(nrow(tmp), replace = T), i]
    }
    return(tmp)
  }
  
  devtools::source_url("https://raw.githubusercontent.com/huayingfang/CCLasso/master/R/SparCC.R")
  
  envir_filter <- enquo(envir_filter)
  
  OTU_ID <- datalist %>%
    with(., if (!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .) %>%
    .$Count_Data %>%
    .$OTU_ID
  
  count <- datalist %>%
    with(., if (!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .) %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    t() %>%
    as_tibble()
  
  if (frac) {
    SparCC_true_cor <- count %>%
      SparCC.frac() %>%
      .$cor.w  
  } else {
    SparCC_true_cor <- count %>%
      SparCC.count() %>%
      .$cor.w  
  }
  
  boot_list <- rep(list(count), n_boot) %>%
    bplapply(., shuffle)
  
  SparCC_boot_cor <- boot_list %>%
    bplapply(., function(x) {
      if (frac) SparCC.frac(x) else SparCC.count(x) }) %>%
    map(., "cor.w") %>%
    map(., function(x) abs(x) >= abs(SparCC_true_cor) & sign(x) == sign(SparCC_true_cor)) %>%
    reduce(`+`)/n_boot
  
  colnames(SparCC_true_cor) <- OTU_ID
  rownames(SparCC_true_cor) <- OTU_ID
  
  return(list(cor = SparCC_true_cor,
              pVal = SparCC_boot_cor))
  
}