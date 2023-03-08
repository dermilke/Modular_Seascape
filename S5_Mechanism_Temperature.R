
library(tidyverse)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

get_mechanism_prop <- function(datalist, bNTI, RC_BC) {
  
  bNTI_mod <- bNTI %>%
    magrittr::set_colnames(rownames(.)) %>%
    as.matrix() 
  
  #bNTI_mod <- bNTI_mod[rownames(RC_BC), rownames(RC_BC)]
  
  bNTI_mod[lower.tri(bNTI_mod)] <- NA
  bNTI_mod[diag(bNTI_mod)] <- NA
  
  RC_BC_mod <- RC_BC %>%
    magrittr::set_rownames(rownames(bNTI_mod)) %>%
    magrittr::set_colnames(colnames(bNTI_mod)) %>%
    as.matrix()
  
  RC_BC_mod[lower.tri(RC_BC_mod)] <- NA
  RC_BC_mod[diag(RC_BC_mod)] <- NA
  
  detailled <- bNTI_mod %>%
    reshape2::melt() %>%
    with(., cbind(., reshape2::melt(RC_BC_mod))) %>%
    magrittr::set_colnames(c("Sample_ID", "To_Sample", "bNTI", "1", "2", "RC_BC")) %>%
    select(Sample_ID, To_Sample, bNTI, RC_BC) %>%
    filter(!is.na(bNTI)) %>%
    mutate(Mechanism = ifelse(bNTI > 2, "Heterogeneous Selection", 
                              ifelse(bNTI < -2, "Homogeneous Selection",
                                     ifelse(RC_BC < 0.95 & RC_BC > -0.95, "Drift",
                                            ifelse(RC_BC > 0.95, "Dispersal Limitation", "Homogenising Dispersal"))))) %>%
    left_join(., select(datalist$Meta_Data, Sample_ID, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("From_Depth_Grp" = "Depth_Grp", "From_Sample" = "Sample_ID", "Sample_ID" = "To_Sample") %>%
    left_join(., select(datalist$Meta_Data, Sample_ID, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("To_Depth_Grp" = "Depth_Grp", "To_Sample" = "Sample_ID") %>%
    filter(From_Depth_Grp == To_Depth_Grp) %>%
    as_tibble()
  
  merged <- detailled %>%
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
  
  return(list(Merged = merged, Detailled = detailled))
  
}

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  #mutate_meta_datalist(Depth_Grp = c(20, 40, 60, 100, 200, 300)[findInterval(Depth, c(20, 40, 60, 100, 200, 300))]) %>%
  mutate_meta_datalist(Ocean = "Atlantic") %>%
  mutate_meta_datalist(Station = as.character(Station)) %>%
  filter_taxa_datalist(Family != "Mitochondria")

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  #mutate_meta_datalist(Depth_Grp = c(20, 40, 60, 100, 200, 300)[findInterval(Depth, c(20, 40, 60, 100, 200, 300))]) %>%
  mutate_meta_datalist(Ocean = "Pacific") %>%
  filter_taxa_datalist(Family != "Mitochondria") %>%
  filter_station_datalist(Depth <= 200)

merged_Atlantic <- bind_rows(
  get_mechanism_prop(datalist_Atlantic, 
                     read.csv("output/Community_Mechanisms/Prokaryotes_weighted_bNTI_FL_Atlantic_complete.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms/Raup_Crick_Atlantic_Prok_FL.csv", row.names = 1))$Detailled %>%
    mutate(Size_Fraction = "0.22-3µm"),
  get_mechanism_prop(datalist_Atlantic, 
                     read.csv("output/Community_Mechanisms/Prokaryotes_weighted_bNTI_SPA_Atlantic_complete.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms/Raup_Crick_Atlantic_Prok_SPA.csv", row.names = 1))$Detailled %>%
    mutate(Size_Fraction = "3-8µm"),
  get_mechanism_prop(datalist_Atlantic, 
                     read.csv("output/Community_Mechanisms/Prokaryotes_weighted_bNTI_LPA_Atlantic_complete.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms/Raup_Crick_Atlantic_Prok_LPA.csv", row.names = 1))$Detailled %>%
    mutate(Size_Fraction = ">8µm")) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c("0.22-3µm", "3-8µm", ">8µm"))) %>%
  #mutate(From_Depth_Grp = ordered(From_Depth_Grp, levels = c(20, 40, 60, 100, 200, 300), labels = c("20m", "40m", "60m", "100m", "200m", "300m"))) %>%
  mutate(Ocean = "Atlantic") 


merged_Pacific <- bind_rows(
  get_mechanism_prop(datalist_Pacific, 
                     read.csv("output/Community_Mechanisms/Prokaryotes_weighted_bNTI_FL_Pacific_complete.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms/Raup_Crick_Pacific_Prok_FL.csv", row.names = 1))$Detailled %>%
    mutate(Size_Fraction = "0.22-3µm"),
  get_mechanism_prop(datalist_Pacific, 
                     read.csv("output/Community_Mechanisms/Prokaryotes_weighted_bNTI_SPA_Pacific_complete.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms/Raup_Crick_Pacific_Prok_SPA.csv", row.names = 1))$Detailled %>%
    mutate(Size_Fraction = "3-8µm"),
  get_mechanism_prop(datalist_Pacific, 
                     read.csv("output/Community_Mechanisms/Prokaryotes_weighted_bNTI_LPA_Pacific_complete.csv", row.names = 1), 
                     read.csv("output/Community_Mechanisms/Raup_Crick_Pacific_Prok_LPA.csv", row.names = 1))$Detailled %>%
    mutate(Size_Fraction = ">8µm")) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c("0.22-3µm", "3-8µm", ">8µm"))) %>%
  #mutate(From_Depth_Grp = ordered(From_Depth_Grp, levels = c(20, 40, 60, 100, 200, 300), labels = c("20m", "40m", "60m", "100m", "200m", "300m"))) %>%
  mutate(Ocean = "Pacific") 

detailled_Atlantic <- merged_Atlantic %>%
  dplyr::rename("From_Sample_ID" = "From_Sample", "To_Sample_ID" = "To_Sample") %>%
  left_join(., datalist_Atlantic$Meta_Data %>% rename_all(function(x) paste0("From_", x)), by = "From_Sample_ID") %>%
  left_join(., datalist_Atlantic$Meta_Data %>% rename_all(function(x) paste0("To_", x)), by = "To_Sample_ID") %>%
  mutate(Temp_Diff = abs(From_Pot_Temperature - To_Pot_Temperature)) %>%
  mutate(Temp_Diff_Grp = seq(0, 28, 4)[findInterval(Temp_Diff, seq(0, 28, 4))]) %>%
  bind_cols(as_tibble(model.matrix(~ Mechanism - 1, .))) %>%
  filter(From_Size_Fraction == To_Size_Fraction) %>%
  filter(From_Depth_Grp.y == To_Depth_Grp.y) %>%
  select(Temp_Diff_Grp, From_Size_Fraction, From_Depth_Grp.y, Mechanism, MechanismDrift, 
         `MechanismDispersal Limitation`, `MechanismHeterogeneous Selection`,
         `MechanismHomogeneous Selection`, `MechanismHomogenising Dispersal`) %>%
  group_by(Temp_Diff_Grp, From_Size_Fraction, From_Depth_Grp.y) %>%
  summarize(MechanismDrift = sum(MechanismDrift),
            `MechanismDispersal Limitation` = sum(`MechanismDispersal Limitation`),
            `MechanismHeterogeneous Selection` = sum(`MechanismHeterogeneous Selection`),
            `MechanismHomogeneous Selection` = sum(`MechanismHomogeneous Selection`),
            `MechanismHomogenising Dispersal` = sum(`MechanismHomogenising Dispersal`)) %>%
  gather("Mechanism", "Value", -Temp_Diff_Grp, -From_Size_Fraction, -From_Depth_Grp.y) 

detailled_Pacific <- merged_Pacific %>%
  dplyr::rename("From_Sample_ID" = "From_Sample", "To_Sample_ID" = "To_Sample") %>%
  left_join(., datalist_Pacific$Meta_Data %>% rename_all(function(x) paste0("From_", x)), by = "From_Sample_ID") %>%
  left_join(., datalist_Pacific$Meta_Data %>% rename_all(function(x) paste0("To_", x)), by = "To_Sample_ID") %>%
  mutate(Temp_Diff = abs(From_Pot_Temperature - To_Pot_Temperature)) %>%
  mutate(Temp_Diff_Grp = seq(0, 28, 4)[findInterval(Temp_Diff, seq(0, 28, 4))]) %>%
  bind_cols(as_tibble(model.matrix(~ Mechanism - 1, .))) %>%
  filter(From_Size_Fraction == To_Size_Fraction) %>%
  filter(From_Depth_Grp.y == To_Depth_Grp.y) %>%
  select(Temp_Diff_Grp, From_Size_Fraction, From_Depth_Grp.y, Mechanism, MechanismDrift, 
         `MechanismDispersal Limitation`, `MechanismHeterogeneous Selection`,
         `MechanismHomogeneous Selection`, `MechanismHomogenising Dispersal`) %>%
  group_by(Temp_Diff_Grp, From_Size_Fraction, From_Depth_Grp.y) %>%
  summarize(MechanismDrift = sum(MechanismDrift),
            `MechanismDispersal Limitation` = sum(`MechanismDispersal Limitation`),
            `MechanismHeterogeneous Selection` = sum(`MechanismHeterogeneous Selection`),
            `MechanismHomogeneous Selection` = sum(`MechanismHomogeneous Selection`),
            `MechanismHomogenising Dispersal` = sum(`MechanismHomogenising Dispersal`)) %>%
  gather("Mechanism", "Value", -Temp_Diff_Grp, -From_Size_Fraction, -From_Depth_Grp.y) 

detailled_Atlantic %>%
  group_by(Temp_Diff_Grp, From_Size_Fraction, From_Depth_Grp.y) %>%
  mutate(Value = Value/sum(Value)) %>%
  mutate(Mechanism = str_remove_all(Mechanism, pattern = "Mechanism")) %>%
  mutate(Mechanism = ordered(Mechanism, levels = c("Homogeneous Selection", "Heterogeneous Selection",
                                                   "Homogenising Dispersal", "Dispersal Limitation", "Drift"),
                             labels = c("Homogenous selection", "Heterogeneous selection",
                                        "Homogenising dispersal", "Dispersal limitation", "Drift"))) %>%
  mutate(From_Size_Fraction = ordered(From_Size_Fraction, levels = c("0.22", "3", "8"),
                                      labels = c("0.22-3 µm", "3-8 µm", "> 8 µm"))) %>%
  mutate(From_Depth_Grp.y = ifelse(From_Depth_Grp.y == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
  
  ggplot(aes(x = Temp_Diff_Grp, y = Value, fill = Mechanism)) +
  geom_bar(stat = "identity", position = "stack", col = "black", size = .2) +
  facet_grid(From_Depth_Grp.y ~ From_Size_Fraction) +
  scale_fill_manual(values = col_palette[c(7,5, 3, 6, 4)]) +
  labs(x = "Temperature difference (°C)", y = "Relative importance of mechanism") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(colour = "black"))

ggsave("figs/Mechanism_Temperature_Atlantic.png", width = 8, height = 5, dpi = 300)


detailled_Pacific %>%
  group_by(Temp_Diff_Grp, From_Size_Fraction, From_Depth_Grp.y) %>%
  mutate(Value = Value/sum(Value)) %>%
  mutate(Mechanism = str_remove_all(Mechanism, pattern = "Mechanism")) %>%
  mutate(Mechanism = ordered(Mechanism, levels = c("Homogeneous Selection", "Heterogeneous Selection",
                                                   "Homogenising Dispersal", "Dispersal Limitation", "Drift"),
                             labels = c("Homogenous selection", "Heterogeneous selection",
                                        "Homogenising dispersal", "Dispersal limitation", "Drift"))) %>%
  mutate(From_Size_Fraction = ordered(From_Size_Fraction, levels = c("0.22", "3", "8"),
                                      labels = c("0.22-3 µm", "3-8 µm", "> 8 µm"))) %>%
  mutate(From_Depth_Grp.y = ifelse(From_Depth_Grp.y == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
  
  ggplot(aes(x = Temp_Diff_Grp, y = Value, fill = Mechanism)) +
  geom_bar(stat = "identity", position = "stack", col = "black", size = .2) +
  facet_grid(From_Depth_Grp.y ~ From_Size_Fraction) +
  scale_fill_manual(values = col_palette[c(7,5, 3, 6, 4)]) +
  labs(x = "Temperature difference (°C)", y = "Relative importance of mechanism") +
  theme_bw() +
  theme(axis.text = element_text(colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(colour = "black"))

ggsave("figs/Mechanism_Temperature_Pacific.png", width = 8, height = 5, dpi = 300)
