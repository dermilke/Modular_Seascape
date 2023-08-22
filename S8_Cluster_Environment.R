#### Cluster Environment Signature #####

source("R/weighted_mean.R")

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_Combined <- combine_data(mutate_meta_datalist(datalist_Atlantic, Station = as.character(Station),
                                                       Ocean = "Atlantic Ocean"), 
                                  mutate_meta_datalist(datalist_Pacific, Station = as.character(Station),
                                                       Ocean = "Pacific Ocean")) %>%
  mutate_meta_datalist(Depth_Int = c(20, 40, 60, 100, 200)[findInterval(Depth, c(20, 40, 60, 100, 200))])

datalist_cluster <- datalist_Combined

datalist_cluster$Count_Data <- datalist_cluster %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., read_csv("output/SparCC/SparCC_Cluster.csv"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  select(-n, -Degree)

cluster_signature <- left_join(weighted_mean_datalist(datalist_cluster, Pot_Temperature),
                               weighted_mean_datalist(datalist_cluster, Salinity),
                               by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, abs(Latitude)),
            by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, Fluorescence),
            by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, Oxygen),
            by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, Depth),
            by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, Si),
            by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, NO3),
            by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, Bacteria_Abundance_Flow),
            by = "Cluster") %>%
  left_join(weighted_mean_datalist(datalist_cluster, Bacterial_Generation_Time),
            by = "Cluster")

p1 <- cluster_signature %>%
  select(-`abs(Latitude)`) %>%
  select_if(is.numeric) %>%
  mutate_all(function(x) x/sum(x)) %>%
  with(., apply(., 1, function(x) (x-mean(x))/sd(x))) %>%
  magrittr::set_rownames(c("Potential Temperature", "Salinity", "Fluorescence",
                           "Oxygen", "Depth", "Silicate", "Nitrate", "Prokaryotic Cell Numbers",
                           "Prokaryotic Generation Time")) %>%
  pheatmap::pheatmap(labels_col = as.character(1:10))

ggsave(p1, filename = "figs/Cluster_Environment.png", width = 6, height = 4, dpi = 300)
