library(tidyverse)

source("R/Import_Data.R")
source("../Jonatan_Project/R/Datalist_Wrangling_Functions.R")
source("R/Import_SparCC_Network.R")
source("R/Similarity_Indices.R")

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

#### Find most important environmental params ####
library(vegan)
library(ggrepel)

params <- c("Pot_Temperature", "Salinity", "Depth", "Fluorescence", "Bacteria_Abundance_Flow", "PO4", "Si", "NO2",
            "HNA_per_LNA", "Eucaryotes_Abundance", "POC", "Bacterial_Generation_Time", "C_N_ratio")

datalist_cluster_mod <- datalist_cluster %>%
  filter_station_datalist(select(datalist_cluster$Meta_Data, params) %>% 
                            with(., rowSums(is.na(.))) == 0)

counts <- datalist_cluster_mod$Count_Data %>%
  select_if(is.numeric) %>%
  t()

model_0 <- cca(counts ~ 1, data = datalist_cluster_mod$Meta_Data %>% select(params))
model_1 <- cca(counts ~ ., data = datalist_cluster_mod$Meta_Data %>% select(params))

ordistep_result <- ordistep(model_0, model_1)

rownames(ordistep_result$CCA$biplot)

new_param_names <- c("Temperature", "Depth", "Salinity", "Prokaryotic cell abundance", "HNA per LNA", "Fluorescence",
                     "POC", "Si", "C to N ratio", "Prokaryotic generation time", "Eukaryotes abundance", "NO2", "PO4")

ggplot(data.frame(ordistep_result$CCA$u), aes(x = CCA1, y = CCA2)) +
  geom_point() +
  geom_segment(aes(x = 0, y = 0, xend = CCA1*2, yend = CCA2*2), col = "darkgrey", size = 1,
               arrow = arrow(length = unit(0.5, "cm")), data = as.tibble(ordistep_result$CCA$biplot) %>% 
                 mutate(Param = rownames(ordistep_result$CCA$biplot))) +
  geom_label_repel(aes(x = CCA1*2, y = CCA2*2, label = Param), data = as.tibble(ordistep_result$CCA$biplot) %>% 
               mutate(Param = new_param_names), size = 3,
               box.padding = .75) +
  geom_point(aes(x = CCA1, y = CCA2, fill = as.factor(Group)), size = 5, data = as.tibble(ordistep_result$CCA$v) %>%
               mutate(Group = 1:10), pch = 21) +
  scale_fill_manual(values = colours) +
  labs(fill = "Module") +
  theme_bw()

ggsave("figs/CCA_Modules.png", width = 6.5, height = 6, dpi = 300)

datalist_Combined_mod <- datalist_Combined %>%
  filter_station_datalist(select(datalist_Combined$Meta_Data, params) %>% 
                            with(., rowSums(is.na(.))) == 0)
  
counts <- datalist_Combined_mod$Count_Data %>%
  select_if(is.numeric) %>%
  t()

model_0 <- cca(counts ~ 1, data = datalist_Combined_mod$Meta_Data %>% select(params))
model_1 <- cca(counts ~ ., data = datalist_Combined_mod$Meta_Data %>% select(params))

ordistep_result_ASVs <- ordistep(model_0, model_1)

rownames(ordistep_result_ASVs$CCA$biplot)

new_param_names <- c("Temperature", "Depth", "HNA per LNA", "Si", "Salinity", "Prokaryotic cell abundance",
                     "PO4", "C to N ratio", "Fluorescence", "POC", "Eukaryotes abundance", "NO2", "Prokaryotic generation time")

module_assignments <- datalist_Combined_mod$Count_Data %>%
  left_join(., read_csv("output/SparCC/SparCC_Cluster.csv"), by = "OTU_ID") %>%
  select(OTU_ID, Cluster, Colour) %>%
  mutate(Colour = ifelse(is.na(Colour), NULL, Colour)) %>%
  mutate(Cluster = ifelse(is.na(Cluster), "Outgroup", Cluster))

Colour_vals <- module_assignments %>%
  select(Cluster, Colour) %>%
  distinct() %>%
  mutate(Cluster = ordered(Cluster, levels = c(1:10, "Outgroup"))) %>%
  arrange(Cluster)

ggplot(data.frame(ordistep_result_ASVs$CCA$v), aes(x = CCA1, y = CCA2)) +
  geom_point(colour = "darkgrey", alpha = .3) +
  geom_point(aes(x = CCA1, y = CCA2, fill = as.factor(Group)), pch = 21, size = 1.8, data = as.tibble(ordistep_result_ASVs$CCA$v) %>%
               mutate(Group = module_assignments$Cluster) %>%
               filter(Group != "Outgroup") %>%
               mutate(Group = ordered(Group, levels = c(1:10)))) +
  scale_fill_manual(values = Colour_vals$Colour) +
  geom_segment(aes(x = 0, y = 0, xend = CCA1*2, yend = CCA2*2), col = "darkgrey", size = 1,
               arrow = arrow(length = unit(0.5, "cm")), data = as.tibble(ordistep_result_ASVs$CCA$biplot)) +
  geom_label_repel(aes(x = CCA1*2, y = CCA2*2, label = Param), data = as.tibble(ordistep_result_ASVs$CCA$biplot) %>% 
                     mutate(Param = new_param_names), size = 3,
                   box.padding = .75) +
  labs(fill = "Module") +
  theme_bw()

ggsave("figs/CCA_ASVs.png", width = 6.5, height = 6, dpi = 300)
