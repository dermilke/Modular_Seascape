library(tidyverse)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

datalist <- import_data("../Malaspina/data/Surface/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) 

datalist_cluster <- datalist

datalist_cluster$Count_Data <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., read_csv("output/SparCC/SparCC_Cluster.csv"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster, Class) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  select(-n, -Degree) %>%
  ungroup()

datatable_cluster <- create_datatable(datalist_cluster, grpBy = Cluster, otherThreshold = 0, 
                                      addOthers = F, addColorScheme = F) 

datatable <- datatable_cluster %>%
  mutate(X_val = ifelse(Variable == "Latitude", Latitude, Longitude)) %>%
  group_by(X_val, Region, Group) %>%
  summarize_if(is.numeric, mean) %>%
  ungroup()

colours <- read_csv("output/SparCC/SparCC_Cluster.csv") %>%
  mutate(Cluster = ordered(Cluster, levels = c(as.character(seq(1, length(unique(Cluster))-1)), "Others"))) %>%
  arrange(Cluster) %>%
  select(Cluster, Colour) %>%
  distinct() %>%
  .$Colour 

p1 <- datatable %>%
  mutate(Group = ordered(Group, levels = c(as.character(1:9), "Others"), labels = c(1:9, 10))) %>%
  filter(Region == "Atlantic_Lat") %>%
  ggplot(., aes(x = X_val, y = Abundance * 100, fill = Group)) +
  geom_area(col = "black", size = .2) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(-25, 35), breaks = seq(-25, 35, by = 15), 
                                    labels = c("25°S", "10°S", "5°N", "20°N", "35°N"), name = "Latitude") +
  labs(fill = "Cluster - Class", y = "Abundance (%)", x = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "none")

p2 <- datatable %>%
  mutate(Group = ordered(Group, levels = c(as.character(1:9), "Others"), labels = c(1:9, 10))) %>%
  filter(Region == "Pacific_Lat") %>%
  ggplot(., aes(x = X_val, y = Abundance * 100, fill = Group)) +
  geom_area(col = "black", size = .2) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(-35, 23), breaks = seq(-35, 25, by = 15), 
                     labels = c("35°S", "20°S", "5°S", "10°N", "25°N"), name = "Latitude") +
  labs(fill = "Cluster - Class", y = "Abundance (%)", x = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "none")

p3 <- datatable %>%
  mutate(Group = ordered(Group, levels = c(as.character(1:9), "Others"), labels = c(1:9, 10))) %>%
  filter(Region == "Transect_North") %>%
  ggplot(., aes(x = X_val, y = Abundance * 100, fill = Group)) +
  geom_area(col = "black", size = .2) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(-155, -20), breaks = seq(-160, -20, by = 40), 
                     labels = c("160°W", "120°W", "80°W", "40°W"), name = "Longitude") +
  labs(fill = "Cluster - Class", y = "Abundance (%)", x = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "none")

p4 <- datatable %>%
  mutate(Group = ordered(Group, levels = c(as.character(1:9), "Others"), labels = c(1:9, 10))) %>%
  filter(Region == "Transect_South") %>%
  ggplot(., aes(x = X_val, y = Abundance * 100, fill = Group)) +
  geom_area(col = "black", size = .2) +
  scale_fill_manual(values = colours) +
  scale_x_continuous(limits = c(-35, 145), breaks = seq(-50, 150, by = 50), 
                     labels = c("50°W", "0°E/W", "50°E", "100°E", "150°E"), name = "Longitude") +
  labs(fill = "Cluster - Class", y = "Abundance (%)", x = "") +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        legend.position = "none")

cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)

ggsave("figs/Network_Cluster_Malaspina.png", width = 9, height = 7, dpi = 300)
