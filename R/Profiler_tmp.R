profiling <- function(datatable, groupFilter) {
  
  tmp <- datatable %>%
    filter(Group == groupFilter)
  
  # 1. Size_Fraction:
  # Only best SF if there are significant differences between SFs
  
  SF_tmp <- group_by(tmp, Size_Fraction) %>% 
    summarize(Abundance = sum(Abundance)) 
  
  if (kruskal.test(tmp$Abundance, tmp$Size_Fraction)$p.value < 0.05) {
    best_SF <- SF_tmp %>%
      arrange(desc(Abundance)) %>%
      dplyr::slice(1) %>%
      .$Size_Fraction
    # Check if best SF is significantly different to other SF
    # and if not: add the other SF that is not different to best SF
    a <- DescTools::DunnTest(tmp$Abundance, as.factor(tmp$Size_Fraction))[[1]][,2] %>%
      tibble::enframe() %>%
      separate(name, c("A", "B"), "-") %>%
      filter(A %in% best_SF | B %in% best_SF) %>%
      filter(value > 0.05)
    
    best_SF <- c(best_SF, a$A[a$A != best_SF], a$B[a$B != best_SF]) %>%
      unique() %>%
      as.character()
    
  } else if (max(SF_tmp$Abundance) / (min(SF_tmp$Abundance)+0.000001) > 10) {
    # Get best SF: If summed Abundance > 10 times smallest: Define best_SF 
    # Then check highest SF_Abundance and look which other SF is at least 33% of Abundance
    # of highest SF Abundance
    best_SF <- SF_tmp %>%
      filter(Abundance > 0) %>%
      filter((Abundance)/max(Abundance) > 0.33) %>%
      select(Size_Fraction) %>%
      as_vector() 
    names(best_SF) <- NULL
    
  } else {
    best_SF <- "All" # or "not 
  }
  
  best_SF <- if (length(best_SF) == 3) "All" else best_SF
  
  # 2. Depths:
  # Best Depth is its 90% Quantile
  best_Depth <- tmp %>%
    {if (best_SF[1] != "All") filter(., Size_Fraction %in% best_SF) else .} %>%
    group_by(Depth) %>%
    summarize(Abundance = sum(Abundance)) %>%
    mutate(Abundance = Abundance / sum(Abundance)) %>%
    with(., quantile(rep(x = Depth, times = round(Abundance*1000)), c(0.1,0.90)))
  
  # 4. Abundance-Maxima:
  
  # First: Calculate a loess-model-fit and find the maxima of the fitted curve. Then take only those maxima
  # that are at least 75% of the highest maximum found. Because it is the fitted curve, outliers wont 
  # affect maximum height very much, thats why i think 15% should be sufficient.
  Model_Fit <- tmp %>%
    {if (best_SF[1] != "All") filter(., Size_Fraction %in% best_SF) else .} %>%
    filter(Depth >= best_Depth[1] & Depth <= best_Depth[2])
  
  if (nrow(Model_Fit) >= 3) {
    
    Model_Fit <- Model_Fit %>%
      mutate(Abundance_Model = loess(Abundance ~ Pot_Temperature, data = ., degree = 2, span = 0.4)$fitted) %>%
      arrange(Pot_Temperature) %>%
      select(Pot_Temperature, Abundance_Model, Abundance) %>%
      distinct() %>%
      mutate(OTU_ID = groupFilter)
    
  }
  
  best_Temp <- tmp %>%
    {if (best_SF[1] != "All") filter(., Size_Fraction %in% best_SF) else .} %>%
    filter(Depth >= best_Depth[1] & Depth <= best_Depth[2]) %>%
    with(., rep(x = Pot_Temperature, times = Abundance*1000)) %>%
    quantile(.,  c(0.1,0.90))
  
  
  
  profile_data <- list(Size_Fraction = best_SF,
                       Depth_Interval = best_Depth,
                       Temp_Interval = best_Temp,
                       Model_Fit = Model_Fit)
  
  return(profile_data)
  
}
