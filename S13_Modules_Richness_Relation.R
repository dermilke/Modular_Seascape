library(tidyverse)

source("R/Import_Data.R")
source("../Jonatan_Project/R/Datalist_Wrangling_Functions.R")
source("R/Import_SparCC_Network.R")
source("R/Similarity_Indices.R")
source("R/Stats_Diversity.R")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_Combined <- combine_data(mutate_meta_datalist(datalist_Atlantic, Station = as.character(Station),
                                                       Ocean = "Atlantic Ocean"), 
                                  mutate_meta_datalist(datalist_Pacific, Station = as.character(Station),
                                                       Ocean = "Pacific Ocean")) %>%
  mutate_meta_datalist(Depth_Int = c(20, 40, 60, 100, 200)[findInterval(Depth, c(20, 40, 60, 100, 200))])

div_table_counts <- datalist_Combined %>%
  rarefy_datalist(rare_lim = 8000, drop = T) %>%
  diversity_datatable()

datalist_cluster <- datalist_Combined %>%
  rarefy_datalist(rare_lim = 8000, drop = T)

datalist_cluster$Count_Data <- datalist_cluster %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., read_csv("output/SparCC/SparCC_Cluster.csv"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  select(-n, -Degree)

div_table_modules <- diversity_datatable(datalist_cluster) 

div_table_modules %>%
  mutate(Richness_Counts = left_join(., div_table_counts, by = "Sample_ID")$Richness.y) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c(0.22, 3, 8), labels = c("0.22-3 µm", "3-8 µm", ">8 µm"))) %>%
  mutate(Depth_Grp = ordered(Depth_Grp, levels = c("Epi", "Meso"), labels = c("20 m - DCM", "DCM - 200 m"))) %>%
  ggplot(aes(x = Richness, y = Richness_Counts, col = Size_Fraction)) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    scale_color_manual(values = cbbPalette[c(4,7,3)]) +
    facet_grid(Ocean~Depth_Grp) +
    labs(y = "ASV richness", x = "No. of modules", col = "Size fraction") +
    theme_bw()

ggsave("figs/Richness_Modules.png", width = 7, height = 6, dpi = 300)

combs <- div_table_modules %>%
  select(Depth_Grp, Size_Fraction, Ocean) %>%
  distinct()

correlation_tests <- pmap(list(x = combs$Depth_Grp, y = combs$Size_Fraction, z = combs$Ocean), function(x, y, z) {
  
  tmp <- div_table_modules %>%
    mutate(Richness_Counts = left_join(., div_table_counts, by = "Sample_ID")$Richness.y) %>%
    filter(Depth_Grp == !!x & Size_Fraction == !!y & Ocean == !!z) %>%
    with(., cor.test(Richness, Richness_Counts))
  
  return(tibble(Depth_Grp = x, Size_Fraction = y, Ocean = z, Cor_Val = tmp$estimate, Cor_P = tmp$p.value))
}) %>%
  bind_rows() %>%
  mutate(Cor_P_adj = p.adjust(Cor_P))


