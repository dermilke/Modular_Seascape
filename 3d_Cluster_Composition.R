library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")

datalist_Reduced <- import_data("data/Combined_Reduced/", kingdom = "Prok", abundance_filter = T, min_counts = 2000)

datalist_cluster <- datalist_Reduced

average_asv_abundance <- datalist_Reduced$Count_Data %>%
  mutate_if(is.numeric, function(x) x/sum(x)*100) %>%
  mutate(average = rowMeans(select_if(., is.numeric))) %>%
  select(OTU_ID, average)

datalist_cluster$Count_Data <- read_csv("output/SparCC/SparCC_Cluster.csv") %>%
  left_join(., average_asv_abundance, by = "OTU_ID") %>%
  select(OTU_ID, Cluster, average) %>%
  reshape2::dcast(OTU_ID~Cluster) %>%
  as_tibble() %>%
  mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
  left_join(., select_if(datalist$Count_Data, is.character), by = "OTU_ID") 
    
list_datatable <- create_datatable(datalist_cluster, grpBy = Family, upper_grp = Class, 
                                   otherThreshold = 0.01, addColorScheme = T)

list_datatable$table %>%
  select(1:4) %>%
  filter(Sample_ID != "Others") %>%
  mutate(Sample_ID = as.numeric(Sample_ID)) %>%
  mutate(Group = plyr::mapvalues(Group, levels(Group), str_replace(levels(Group), pattern = ";", replacement = " - "))) %>%
  mutate(Group = ordered(Group, levels = levels(datatable$Group), 
                         labels = str_replace_all(levels(datatable$Group), pattern = "Unknown Marinimicrobia \\(SAR406 clade\\) \\- Unknown Marinimicrobia \\(SAR406 clade\\)",
                                                  replacement = "Unknown Marinimicrobia (SAR406 clade)"))) %>%

  ggplot(., aes(x = as.factor(Sample_ID), y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = "stack", width = .8, color = "black", size = .3) +
    scale_fill_manual(values = list_datatable$color) +
    labs(fill = "Class - Family", x = "Cluster", y = "Relative abundance within dataset (%)") +
    theme_bw() +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          legend.key.size = unit(.8, "line")) +
    guides(fill = guide_legend(nrow = ceiling(length(unique(datatable$Group))), 
                               title.position = "top"))

ggsave("figs/Network_Cluster_Composition_rel_dataset.png", width = 10, height = 5.5, dpi = 300)
