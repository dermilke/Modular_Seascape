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

network_properties <- read_csv("output/Network_Properties.csv") %>%
  filter(Type == "Positive")

div_table_modules <- diversity_datatable(datalist_cluster) %>%
  mutate(Modularity = left_join(datalist_cluster$Meta_Data, 
                                network_properties %>% select(Sample_ID, Modularity), by = "Sample_ID")$Modularity)
  
Div_table_tmp <- div_table_modules %>%
  mutate(Ocean = ifelse(Ocean == "Pacific Ocean", "Pacific", "Atlantic")) %>%
  mutate(Latitude = round(Latitude, digits = 2)) %>%
  select(Station, Latitude, Sample_ID, Ocean, Depth, Size_Fraction, Depth_Grp, Richness, Modularity, ESN) %>%
  distinct()
  
tracer_concentration <- full_join(Div_table_tmp, 
                                  read_csv("output/Tracer_Concentration.csv") %>%
                                    mutate(Latitude = round(Latitude, digits = 2)), by = c("Latitude", "Ocean")) %>%
  mutate(Tracer_Concentration = ifelse(is.nan(Tracer_Concentration), 0, Tracer_Concentration)) %>%
  filter(Depth_Grp == "Epi")

bind_rows(select(tracer_concentration, Latitude, Tracer_Concentration, Ocean, Year_Simulation, Size_Fraction, Depth_Grp) %>%
            mutate(Type = "Tracer concentration") %>%
            dplyr::rename("Number" = "Tracer_Concentration"),
          select(tracer_concentration, Latitude, ESN, Ocean, Year_Simulation, Size_Fraction, Depth_Grp) %>%
            mutate(Type = "Effective number\nof modules") %>%
            dplyr::rename("Number" = "ESN")) %>%
  bind_rows(., select(tracer_concentration, Latitude, Modularity, Ocean, Year_Simulation, Size_Fraction, Depth_Grp) %>%
              mutate(Type = "Modularity") %>%
              dplyr::rename("Number" = "Modularity")) %>%
  filter(Year_Simulation == 5) %>%
  filter((Type == "Tracer concentration" & Size_Fraction == 0.22) | (Type != "Tracer concentration")) %>%
  filter(Latitude > -60) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c(0.22, 3, 8), labels = c("0.22-3 µm", "3-8 µm", ">8 µm"))) %>%
  mutate(Type = ifelse(Type == "Tracer concentration", "Tracer concentration\n(simulation time: 5 years)", Type)) %>%
  ggplot(aes(x = Latitude, y = Number, col = Type, shape = as.factor(Size_Fraction))) +
    geom_point() +
    geom_smooth(span = .35, se = F) +
    facet_grid(Type~Ocean, scales = "free") +
    scale_x_continuous(breaks = c(-60, -30, 0, 30, 60), 
                       labels = c("60°S", "30°S", "0°S/N", "30°N", "60°N"), limits = c(-60, 60)) +
    scale_color_manual(values = cbbPalette[c(3, 2, 4)]) +
    labs(col = "", shape = "Size fraction", y = "") +
    theme_bw()

ggsave("figs/Tracer_concentration_network_year_5.png", width = 9, height = 8, dpi = 300)
