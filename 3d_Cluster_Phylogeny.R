#### Phylogenetic relatedness of Cluster ####
library(ggtree)
library(doParallel)

datalist <- import_data("data/Combined/", kingdom = "Prok", abundance_filter = T, min_counts = 2000)

datalist_cluster <- datalist

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

Tree <- ape::read.tree("output/FastTree_Tree/Prok_Combined_SparCC.tree")

Count_Table <- datalist_cluster$Count_Data %>%
  select_if(is.numeric) %>%
  as.matrix() %>%
  magrittr::set_rownames(datalist_cluster$Count_Data$OTU_ID) %>%
  Matrix::Matrix()

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

ggsave("figs/Phylo_Dist_Modules_Tree.png", width = 5, height = 5, dpi = 300)