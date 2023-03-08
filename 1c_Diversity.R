library(tidyverse)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")
source("R/Stats_Diversity.R")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Atlantic") %>%
  mutate_meta_datalist(Station = as.character(Station)) %>%
  filter_taxa_datalist(Family != "Mitochondria")

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Pacific") %>%
  filter_taxa_datalist(Family != "Mitochondria")

datalist_Combined <- combine_data(datalist_Atlantic, datalist_Pacific)

diversity <- datalist_Combined %>%
  mutate_meta_datalist(Sigma = gsw::gsw_sigma0(Salinity, gsw::gsw_CT_from_pt(Salinity, Pot_Temperature))) %>%
  rarefy_datalist(rare_lim = 8000, drop = T) %>%
  diversity_datatable() %>%
  arrange(Ocean)

diversity %>%
  filter(Depth_Grp == "Epi") %>%
  mutate(Size_Fraction = ordered(Size_Fraction, labels = c("0.22-3 µm", "3-8 µm", ">8 µm"))) %>%
  mutate(Depth_Grp = ifelse(Depth_Grp == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
  
  ggplot(., aes(x = Latitude, y = Richness, col = Size_Fraction)) +
    geom_point(size = 2) +
    geom_smooth(aes(col = Size_Fraction), se = F, size = 1, span = .4) +
    scale_y_log10() +
    scale_color_manual(values = cbbPalette[c(4,7,3)]) +
    labs(y = "Richness", col = "Size fraction") +
    facet_grid(Ocean~Depth_Grp) +
    scale_x_continuous(limits = c(-63, 60), breaks = seq(-60, 60, 30), 
                       labels = c("60°S", "30°S", "0°S/N", "30°N", "60°N")) +
    theme_bw()

ggsave("figs/Diversity/Richness_Epi.png", width = 5, height = 6, dpi = 300)