import_sparcc_network <- function(cor_file, pval_file, min_r, min_p) {
  
  cor <- read.csv(cor_file, header = T,
                  check.names = F) %>%
    magrittr::set_rownames(.[,1]) %>%
    .[,-1] %>%
    as.matrix()
  
  pval <- read.csv(pval_file, check.names = F) %>%
    magrittr::set_rownames(.[,1]) %>%
    .[,-1] %>%
    as.matrix() %>%
    p.adjust() %>%
    matrix(., nrow = nrow(cor), ncol = ncol(cor), byrow = T)
  
  cor[upper.tri(cor)] <- 0
  cor[abs(cor) < min_r] <- 0
  cor[pval > min_p] <- 0
  
  edges <- cor %>%
    reshape2::melt(as.is = T) %>%
    dplyr::rename("From" = "Var1", "To" = "Var2", "weight" = "value") %>%
    as_tibble() %>%
    filter(weight != 0) %>%
    mutate(Type = ifelse(weight < 0, "Negative", "Positive"))
  
  return(edges)
  
}