library(igraph)
library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Atlantic") 

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Pacific") 

datalist <- combine_data(mutate_meta_datalist(datalist_Atlantic, Station = as.character(Station)),
                         mutate_meta_datalist(datalist_Pacific, Station = as.character(Station)))

network_pos <- read_graph("output/SparCC_Network_pos.txt", format = "graphml")

results <- tibble(Sample_ID = NA, Modularity = NA, Transitivity = NA, Mean_Dist = NA, Type = NA)

for (i in datalist$Meta_Data$Sample_ID) {
  
  tmp <- datalist %>%
    filter_station_datalist(Sample_ID == !!i)

  sub_network_pos <- subgraph(network_pos, which(names(V(network_pos)) %in% tmp$Count_Data$OTU_ID)) %>%
    delete_vertices(., which(degree(.) == 0))
  
  tmp_modul_pos <- cluster_walktrap(sub_network_pos, weights = NULL) %>%
    modularity()
  
  tmp_transit_pos <- transitivity(sub_network_pos, type = "average")
  
  tmp_dist_pos <- mean_distance(sub_network_pos)
  
  results <- results %>% 
    add_case(Sample_ID = !!i, Modularity = tmp_modul_pos, Transitivity = tmp_transit_pos, 
             Mean_Dist = tmp_dist_pos, Type = "Positive") 
}

data <- results %>%
  dplyr::slice(-1) %>%
  left_join(., datalist$Meta_Data, by = "Sample_ID")

write_csv(data, "output/Network_Properties.csv")

data <- read_csv("output/Network_Properties.csv")

p1 <- data %>%
  filter(Type == "Positive") %>%
  mutate(Kombi = paste0(Ocean, " - ", Size_Fraction)) %>%
  mutate(Depth_Grp = ifelse(Depth_Grp == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
ggplot(., aes(x = Latitude, y = Modularity, col = as.factor(Size_Fraction))) +
  geom_point() +
  geom_smooth(span = .35, se = F) +
  scale_color_manual(values = cbbPalette[c(4,7,3)]) +
  scale_x_continuous(breaks = c(-60, -30, 0, 30, 60), 
                     labels = c("60°S", "30°S", "0°S/N", "30°N", "60°N"), limits = c(-60, 60)) +
  facet_grid(Depth_Grp ~ Ocean) +
  labs(col = "Size fraction")

legend <- cowplot::get_legend(p1) 

p1 <- p1 +
  theme_bw() +
  theme(legend.position = "none")

p2 <- data %>%
  filter(Type == "Positive") %>%
  mutate(Kombi = paste0(Ocean, " - ", Size_Fraction)) %>%
  mutate(Depth_Grp = ifelse(Depth_Grp == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
  filter(Transitivity > 0.4) %>%
  ggplot(., aes(x = Latitude, y = Transitivity, col = as.factor(Size_Fraction))) +
  geom_point() +
  geom_smooth(span = .35, se = F) +
  scale_color_manual(values = cbbPalette[c(4,7,3)]) +
  scale_x_continuous(breaks = c(-60, -30, 0, 30, 60), 
                     labels = c("60°S", "30°S", "0°S/N", "30°N", "60°N"), limits = c(-60, 60)) +
  facet_grid(Depth_Grp ~ Ocean) +
  labs(col = "Size fraction", y = "Clustering coefficient") +
  theme_bw() +
  theme(legend.position = "none")

cowplot::plot_grid(p2, p1, legend, nrow = 1, ncol = 3)

ggsave("figs/network_properties_pos.pdf", width = 16, height = 5, dpi = 300)
