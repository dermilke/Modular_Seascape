library(igraph)
library(tidyverse)

source("R/Import_Data.R")
source("../Jonatan_Project/R/Datalist_Wrangling_Functions.R")
source("R/Import_SparCC_Network.R")

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist <- combine_data(mutate_meta_datalist(datalist_Atlantic, Station = as.character(Station),
                                              Ocean = "Atlantic Ocean"), 
                         mutate_meta_datalist(datalist_Pacific, Station = as.character(Station),
                                              Ocean = "Pacific Ocean"))

SparCC_Cluster <- read_csv("output/SparCC/SparCC_Cluster.csv")

sequences <- seqinr::read.fasta("data/Combined_Reduced/Fasta/Prok/Prok_Subset_SparCC.fasta", as.string = T,
                                forceDNAtolower = F, set.attributes = F)

enframe(sequences[SparCC_Cluster$OTU_ID]) %>%
  mutate(Sequence = unlist(value)) %>%
  select(-value) %>%
  dplyr::rename("OTU_ID" = "name") %>%
  left_join(., datalist$Count_Data %>% select_if(is.character), by = "OTU_ID") %>%
  left_join(., select(SparCC_Cluster, OTU_ID, Cluster), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  dplyr::rename("ASV_hash" = "OTU_ID") %>%
  filter(Family != "Mitochondria") %>%
  filter(!is.na(Kingdom)) %>%
  write_csv("output/Network_Modules_Table.csv")
