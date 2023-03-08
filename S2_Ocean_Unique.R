library(tidyverse)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

get_colors <- function(datatable) {
  colorTaxonomy <- read_csv("~/PhD/SoftwareBuilds/ExCom/Data/Colors/Taxonomy_Colour.csv")
  
  derivative_color = str_split_fixed(datatable$Group, " - ", 2) %>%
    unique() %>%
    .[,1] %>%
    table() %>%
    as_tibble() %>%
    dplyr::rename("Group" = ".") %>%
    mutate(Group = ordered(Group, levels = c(unique(Group)[-which(Group == "Others")], "Others"))) %>%
    arrange(Group) %>%
    mutate(main_color = colorTaxonomy$Colour[match(Group, colorTaxonomy$Group)]) %>%
    mutate(rowNumber = 1:n()) %>%
    with(., {
      derivative_color <- NULL
      for (i in rowNumber) {
        
        derivative_color <- c(derivative_color,
                              colorspace::lighten(main_color[i], 
                                                  amount = {if (n[i] == 1) 0
                                                    else if (n[i] > 6) (seq_len(n[i])/n[i])*1.6-0.8
                                                    else (seq_len(n[i])/n[i])-0.5}
                              )
        )
        
      }
      derivative_color[length(derivative_color)] <- "grey40"
      derivative_color
    })
  derivative_color
}

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Atlantic") %>%
  mutate_meta_datalist(Station = as.character(Station)) %>%
  filter_taxa_datalist(Family != "Mitochondria") %>%
  filter_station_datalist(Depth <= 200)

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Pacific") %>%
  filter_taxa_datalist(Family != "Mitochondria") %>%
  filter_station_datalist(Depth <= 200)

datalist_Combined <- combine_data(datalist_Atlantic, datalist_Pacific)

datalist_summarized <- datalist_Combined %>%
  summarize_by_param(Ocean, Latitude, Size_Fraction, Depth_Grp)

datalist_mod <- datalist_summarized 

ocean_unique <- datalist_mod$Count_Data %>%
  select_if(is.numeric) %>%
  mutate_if(datalist_mod$Meta_Data$Ocean == "Atlantic", function(x) ifelse(x > 0, 1, NA)) %>%
  mutate_if(datalist_mod$Meta_Data$Ocean == "Pacific", function(x) ifelse(x > 0, 2, NA)) %>%
  with(., apply(., 1, mean, na.rm = T))

datalist_mod$Count_Data <- datalist_mod$Count_Data %>%
  mutate(Ocean_Taxa = ifelse(ocean_unique == 1, "Atlantic", ifelse(ocean_unique == 2, "Pacific", "Mix"))) 
  
datatable <- datalist_mod %>%
  filter_taxa_datalist(Ocean_Taxa != "Mix") %>%
  create_datatable(grpBy = Family, upper_grp = Class, otherThreshold = 0.01,
                   addColorScheme = T, addOthers = T) 

datatable$table %>%
  mutate(Depth_Grp = ifelse(Depth_Grp == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c(0.22, 3, 8), 
                                 labels = c("0.22-3 µm", "3-8 µm", ">8 µm"))) %>%
  mutate(Kombi = paste0(Ocean, "\n", Depth_Grp)) %>%
ggplot(., aes(x = as.factor(round(Latitude, digits = 2)), y = Abundance * 100, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = datatable$color) +
  facet_grid(Size_Fraction~Kombi, scales = "free_x") +
  labs(y = "Abundance (%)", x = "Latitude", fill = "Class - Family") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.text = element_text(size = 5)) +
  guides(fill = guide_legend(ncol = 1, keywidth = .5, keyheight = .5))

ggsave(filename = "figs/Ocean_Unique_ASVs_final.png", width = 9, height = 6, dpi = 300)