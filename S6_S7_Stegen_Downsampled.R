# Get relative importance of ecological mechanisms from Stegen et al. framework
# Requires bNTI table and RC_BC matrix
# -> See R-scripts folder

library(tidyverse)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

get_mechanism_prop <- function(datalist, bNTI, RC_BC) {
  
  datalist_tmp <- datalist %>%
    filter_station_datalist(Sample_ID %in% rownames(!!bNTI) &
                            Sample_ID %in% rownames(!!RC_BC))
  
  bNTI_mod <- bNTI %>%
    magrittr::set_colnames(rownames(.)) %>%
    as.matrix() %>%
    .[datalist_tmp$Meta_Data$Sample_ID, 
      datalist_tmp$Meta_Data$Sample_ID]
  
  bNTI_mod[lower.tri(bNTI_mod)] <- NA
  bNTI_mod[diag(bNTI_mod)] <- NA
  
  RC_BC_mod <- RC_BC %>%
    magrittr::set_colnames(rownames(RC_BC)) %>%
    as.matrix() %>%
    .[datalist_tmp$Meta_Data$Sample_ID, 
      datalist_tmp$Meta_Data$Sample_ID]
  
  RC_BC_mod[lower.tri(RC_BC_mod)] <- NA
  RC_BC_mod[diag(RC_BC_mod)] <- NA
  
  merged <- bNTI_mod %>%
    reshape2::melt() %>%
    with(., cbind(., reshape2::melt(RC_BC_mod))) %>%
    magrittr::set_colnames(c("Sample_ID", "To_Sample", "bNTI", "1", "2", "RC_BC")) %>%
    select(Sample_ID, To_Sample, bNTI, RC_BC) %>%
    filter(!is.na(bNTI)) %>%
    mutate(Mechanism = ifelse(bNTI > 2, "Heterogeneous Selection", 
                              ifelse(bNTI < -2, "Homogeneous Selection",
                                     ifelse(RC_BC < 0.95 & RC_BC > -0.95, "Drift",
                                            ifelse(RC_BC > 0.95, "Dispersal Limitation", "Homogenising Dispersal"))))) %>%
    left_join(., select(datalist_tmp$Meta_Data, Sample_ID, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("From_Depth_Grp" = "Depth_Grp", "From_Sample" = "Sample_ID", "Sample_ID" = "To_Sample") %>%
    left_join(., select(datalist_tmp$Meta_Data, Sample_ID, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("To_Depth_Grp" = "Depth_Grp", "To_Sample" = "Sample_ID") %>%
    filter(From_Depth_Grp == To_Depth_Grp) %>%
    group_by(Mechanism, From_Depth_Grp) %>%
    summarize(Num = n()) %>%
    group_by(From_Depth_Grp) %>%
    mutate(Num = Num/sum(Num)) %>%
    ungroup() %>%
    mutate(Mechanism_Group = ordered(ifelse(Mechanism == "Homogeneous Selection" | Mechanism == "Heterogeneous Selection", "Selection",
                                            ifelse(Mechanism == "Homogenising Dispersal" | Mechanism == "Dispersal Limitation", "Dispersal",
                                                   "Drift")),
                                     levels = c("Drift", "Dispersal", "Selection"))) %>%
    mutate(Mechanism = ordered(Mechanism, levels = c("Homogeneous Selection", "Heterogeneous Selection", 
                                                     "Homogenising Dispersal", "Dispersal Limitation", "Drift")))
  
  return(merged)
  
}

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Atlantic") %>%
  mutate_meta_datalist(Station = as.character(Station)) %>%
  filter_taxa_datalist(Family != "Mitochondria")

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Pacific") %>%
  filter_taxa_datalist(Family != "Mitochondria") %>%
  filter_station_datalist(Depth <= 200)

datalist_Combined <- import_data("data/Combined_Reduced/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = ifelse(str_detect(Cruise, pattern = "^ANT"), "Atlantic", "Pacific")) %>%
  filter_taxa_datalist(Family != "Mitochondria") %>%
  filter_station_datalist(Depth <= 200) #%>%
  mutate_meta_datalist(Depth_Grp = c(20, 40, 60, 100, 200)[findInterval(Depth, c(20, 40, 60, 100, 200))])

datalist_Combined_reduced <- datalist_Combined %>%
  filter_station_datalist(!(Station %in% c("17_2_P2", "20_P2", "297", "8_2_P2", "6_P2", "307", "310", "313", "318", "322", "326")))

datalist_Atlantic_reduced <- datalist_Combined_reduced %>%
  filter_station_datalist(Ocean == "Atlantic")

datalist_Pacific_reduced <- datalist_Combined_reduced %>%
  filter_station_datalist(Ocean == "Pacific")

merged_Atlantic <- bind_rows(
  get_mechanism_prop(datalist = datalist_Atlantic_reduced,
                     bNTI = read.csv("output/Community_Mechanisms_Downsampled/Prokaryotes_weighted_bNTI_FL_Atlantic_downsampled.csv", row.names = 1), 
                     RC_BC = read.csv("output/Community_Mechanisms_Downsampled/Raup_Crick_Prok_FL_Atlantic_Downsampled.csv", row.names = 1)) %>%
    mutate(Size_Fraction = "0.22-3µm"),
  get_mechanism_prop(datalist_Atlantic_reduced, 
                     read.csv("output/Community_Mechanisms_Downsampled/Prokaryotes_weighted_bNTI_SPA_Atlantic_downsampled.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms_Downsampled/Raup_Crick_Prok_SPA_Atlantic_Downsampled.csv", row.names = 1)) %>%
    mutate(Size_Fraction = "3-8µm"),
  get_mechanism_prop(datalist_Atlantic_reduced, 
                     read.csv("output/Community_Mechanisms_Downsampled/Prokaryotes_weighted_bNTI_LPA_Atlantic_downsampled.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms_Downsampled/Raup_Crick_Prok_LPA_Atlantic_Downsampled.csv", row.names = 1)) %>%
    mutate(Size_Fraction = ">8µm")) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c("0.22-3µm", "3-8µm", ">8µm"))) %>%
  mutate(Ocean = "Atlantic")

merged_Pacific <- bind_rows(
  get_mechanism_prop(datalist = datalist_Pacific_reduced,
                     bNTI = read.csv("output/Community_Mechanisms_Downsampled/Prokaryotes_weighted_bNTI_FL_Pacific_downsampled.csv", row.names = 1), 
                     RC_BC = read.csv("output/Community_Mechanisms_Downsampled/Raup_Crick_Prok_FL_Pacific_Downsampled.csv", row.names = 1)) %>%
    mutate(Size_Fraction = "0.22-3µm"),
  get_mechanism_prop(datalist_Pacific_reduced, 
                     read.csv("output/Old_Community_Mechanisms/Prokaryotes_weighted_bNTI_SPA_FastTree.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms_Downsampled/Raup_Crick_Prok_SPA_Pacific_Downsampled.csv", row.names = 1)) %>%
    mutate(Size_Fraction = "3-8µm"),
  get_mechanism_prop(datalist_Pacific_reduced, 
                     read.csv("output/Old_Community_Mechanisms/Prokaryotes_weighted_bNTI_LPA_FastTree.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms_Downsampled/Raup_Crick_Prok_LPA_Pacific_Downsampled.csv", row.names = 1)) %>%
    mutate(Size_Fraction = ">8µm")) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c("0.22-3µm", "3-8µm", ">8µm"))) %>%
  mutate(Ocean = "Pacific") 

bind_rows(merged_Atlantic, merged_Pacific) %>%
  mutate(Kombi = paste0(Ocean, "\n", Size_Fraction)) %>%
  mutate(Kombi = ordered(Kombi, levels = c("Atlantic\n0.22-3µm", "Atlantic\n3-8µm", "Atlantic\n>8µm",
                                           "Pacific\n0.22-3µm", "Pacific\n3-8µm", "Pacific\n>8µm"))) %>%
  #filter(From_Depth_Grp == "Epi") %>%
  #mutate(From_Depth_Grp = ifelse(From_Depth_Grp == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
  
  ggplot(., aes(x = Mechanism_Group, y = Num*100, fill = Mechanism)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(y = "Proportion of Mechanisms (%)", x = "") +
  scale_fill_manual(values = cbbPalette[c(7,5, 3, 6, 4)]) +
  facet_grid(Kombi~From_Depth_Grp) +
  theme_bw() +
  theme(legend.position="right") +
  guides(fill = guide_legend(nrow = 5, title.position = "top")) +
  ylim(c(0, 75))

ggsave("figs/Shaping_Mechanisms_Comparison_Downsampled_Groups.png", width = 7, height = 8, dpi = 300)
