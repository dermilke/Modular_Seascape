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

n <- read_csv("output/SparCC/SparCC_Cluster.csv") %>%
  select(Cluster, n) %>%
  distinct() %>%
  group_by(Cluster) %>%
  summarize_all(sum) %>%
  arrange(Cluster) %>%
  .$n

datatable_cluster <- create_datatable(datalist_cluster, grpBy = Cluster, upper_grp = Cluster,
                                      otherThreshold = 0) %>%
  group_by(Ocean, Latitude, Depth_Grp, Size_Fraction, Group) %>%
  summarize_if(is.numeric, mean) %>%
  ungroup() %>%
  group_by(Group)

colours <- read_csv("output/SparCC/SparCC_Cluster.csv") %>%
  mutate(Cluster = ordered(Cluster, levels = c(as.character(seq(1, length(unique(Cluster))-1)), "Others"))) %>%
  arrange(Cluster) %>%
  select(Cluster, Colour) %>%
  distinct() %>%
  .$Colour 

# --> Change size-fraction to one of c(0.22, 3, 8) for different plots

datatable_cluster %>%
  mutate(Depth_Grp = ifelse(Depth_Grp == "Epi", "20 m - DCM", "DCM - 200 m")) %>%
  mutate(Group = ordered(Group, levels = c(as.character(1:13), "Others"))) %>%
  filter(Size_Fraction == 0.22) %>%
  
  ggplot(., aes(x = Latitude, y = Abundance * 100, fill = Group)) +
    geom_area(col = "black", size = .2) +
    scale_fill_manual(values = colours) +
    facet_grid(Depth_Grp~Ocean) +
    scale_x_continuous(limits = c(-63, 60), breaks = seq(-60, 60, 30),
                       labels = c("60°S", "30°S", "0°S/N", "30°N", "60°N")) +
    labs(fill = "Cluster", y = "Abundance (%)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

ggsave("figs/Network_Cluster_Abundance_FL_Complete.png", width = 8, height = 7, dpi = 300)

#### Get taxonomic composition of clusters with average abundances in dataset ####

average_asv_abundance <- datalist$Count_Data %>%
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

datatable <- list_datatable$table %>%
  select(1:4) %>%
  filter(Sample_ID != "Others") %>%
  mutate(Sample_ID = as.numeric(Sample_ID)) %>%
  mutate(Group = plyr::mapvalues(Group, levels(Group), str_replace(levels(Group), pattern = ";", replacement = " - ")))

colorvalues <- list_datatable$color

Group_Levels <- str_replace_all(levels(datatable$Group), pattern = "Unknown Marinimicrobia \\(SAR406 clade\\) \\- Unknown Marinimicrobia \\(SAR406 clade\\)",
                replacement = "Unknown Marinimicrobia (SAR406 clade)")

datatable %>%
  mutate(Group = ordered(Group, levels = levels(datatable$Group), labels = Group_Levels)) %>%
  
  ggplot(., aes(x = as.factor(Sample_ID), y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = "stack", width = .8, color = "black", size = .3) +
    scale_fill_manual(values = colorvalues) +
    labs(fill = "Class - Family", x = "Cluster", y = "Relative abundance within dataset (%)") +
    theme_bw() +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          legend.key.size = unit(.8, "line")) +
    guides(fill = guide_legend(nrow = ceiling(length(unique(datatable$Group))), 
                               title.position = "top"))

ggsave("figs/Network_Cluster_Composition_rel_dataset.png", width = 10, height = 5.5, dpi = 300)

#### Phylogenetic relatedness of Cluster ####
library(ggtree)

Tree <- ape::read.tree("output/FastTree_Tree/Prok_Combined_SparCC.tree")

Count_Table <- datalist_cluster$Count_Data %>%
  select_if(is.numeric) %>%
  as.matrix() %>%
  Matrix::Matrix() %>%
  magrittr::set_rownames(datalist_cluster$Count_Data$OTU_ID)

my.ps <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(Count_Table), taxa_are_rows=T), 
                            phyloseq::phy_tree(Tree))

Unifrac_dist <- distance_wrapper(my.ps, method = "UniFrac_weighted", size.thresh = 1, pseudocount = 10^-6, nblocks = 100, 
                                 use.cores = 4, cor.use = "na.or.complete")

Unifrac_dist %>%
  magrittr::set_colnames(c(1:10)) %>%
  magrittr::set_rownames(c(1:10)) %>%
  as.dist() %>%
  hclust(method = "mcquitty") %>%
  plot(., hang = -1, main = "", xlab = "", ylab = "", axes = F)
