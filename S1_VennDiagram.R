#### VennDiagrams ####

library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

datalist <- import_data("data/Combined_Reduced/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_summarized <- datalist %>%
  summarize_by_param(Ocean)
  
datalist_summarized_Depth_Grp <- datalist %>%
  summarize_by_param(Depth_Grp)

datalist_summarized_SF <- datalist %>%
  summarize_by_param(Size_Fraction)

Depth_Grp <- datalist_summarized_Depth_Grp$Count_Data %>%
  select_if(is.numeric) %>%
  mutate_all(function(x) x>0) %>%
  mutate(Depth_Grp = ifelse(rowSums(.) == 2, "Both depthlayer", 
                            ifelse(Epi == TRUE, "Epipelagic", "Mesopelagic"))) %>%
  .$Depth_Grp

Size_Fraction <- datalist_summarized_SF$Count_Data %>%
  select_if(is.numeric) %>%
  mutate_all(function(x) x>0) %>%
  magrittr::set_names(c("FL", "SPA", "LPA"))

logic_mat <- datalist_summarized$Count_Data %>%
  select_if(is.numeric) %>%
  mutate_all(function(x) x>0) %>%
  mutate(Depth_Grp = Depth_Grp) %>%
  bind_cols(., Size_Fraction)

venn_data_022 <- Venn(list(
  Pacific = datalist$Count_Data$OTU_ID[logic_mat$Pacific & logic_mat$FL], 
  Atlantic = datalist$Count_Data$OTU_ID[logic_mat$Atlantic & logic_mat$FL]
  )) %>%
  process_data()

venn_data_3 <- Venn(list(
  Pacific = datalist$Count_Data$OTU_ID[logic_mat$Pacific & logic_mat$SPA], 
  Atlantic = datalist$Count_Data$OTU_ID[logic_mat$Atlantic & logic_mat$SPA]
)) %>%
  process_data()

venn_data_8 <- Venn(list(
  Pacific = datalist$Count_Data$OTU_ID[logic_mat$Pacific & logic_mat$LPA], 
  Atlantic = datalist$Count_Data$OTU_ID[logic_mat$Atlantic & logic_mat$LPA]
)) %>%
  process_data()
  
venn_plotter <- function(venn_obj) {
  ggplot() +
    geom_sf(aes(fill = name), data = venn_region(venn_obj)) +
    geom_sf(size = .5, color = "black", data = venn_setedge(venn_obj), show.legend = F) +
    geom_sf_text(aes(label = name), data = venn_setlabel(venn_obj), fontface = "bold", size = 5) +
    geom_sf_label(aes(label=paste0(round(count/sum(count)*100, digits = 1), "%\n", count)), 
                  fontface = "bold", data = venn_region(venn_obj),
                  col = "black", size = 5) +
    scale_fill_manual(values = cbbPalette[c(3, 3, 4)]) +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA))
}

p1 <- venn_plotter(venn_data_022) +
  theme(axis.title.y = element_text(size = 20, angle = 90, face = "bold")) +
  labs(y = "0.22-3 µm")
p2 <- venn_plotter(venn_data_3) +
  theme(axis.title.y = element_text(size = 20, angle = 90, face = "bold")) +
  labs(y = "3-8 µm")
p3 <- venn_plotter(venn_data_8) +
  theme(axis.title.y = element_text(size = 20, angle = 90, face = "bold")) +
  labs(y = ">8 µm")

cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave("figs/Venn_Diagram_Oceans.png", dpi = 300, width = 4, height = 8)  


