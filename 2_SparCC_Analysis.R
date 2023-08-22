library(igraph)
library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")
source("R/Import_SparCC_Network.R")

get_colors_cont <- function(vec, palette, reverse = F, n = 9) {
  
  vec <- ifelse(is.na(vec), 0, vec)
  
  if (reverse) { ramp <- rev(RColorBrewer::brewer.pal(n, palette)) %>% colorRamp(.) } else 
  { ramp <- RColorBrewer::brewer.pal(n, palette) %>%  colorRamp(.) } 
  
  ramp(vegan::decostand(vec, method = "range")) %>%
    rgb(., maxColorValue = 255)
}

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso"))

datalist_Combined <- combine_data(mutate_meta_datalist(datalist_Atlantic, Station = as.character(Station),
                                                       Ocean = "Atlantic Ocean"), 
                                  mutate_meta_datalist(datalist_Pacific, Station = as.character(Station),
                                                       Ocean = "Pacific Ocean")) %>%
  mutate_meta_datalist(Depth_Int = c(20, 40, 60, 100, 200)[findInterval(Depth, c(20, 40, 60, 100, 200))])

Max_Count <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  mutate(Max = apply(., 1, which.max)) %>%
  mutate(OTU_ID = datalist$Count_Data$OTU_ID) %>%
  select(OTU_ID, Max) %>%
  cbind(., datalist$Meta_Data[.$Max,]) %>%
  as_tibble() %>%
  left_join(., select_if(datalist$Count_Data, is.character), by = "OTU_ID")

r_threshold = 0.51

files <- list.files("Output/SparCC", pattern = "^Cor_SparCC_Prok.*", full.name = T)

data_combined <- map(files, function(x) {
  import_sparcc_network(cor_file = x, 
                        pval_file = str_replace_all(x, pattern = "Cor", replacement = "Pval"),
                        min_r = r_threshold, min_p = 0.05)
}) %>%
  bind_rows() %>%
  filter(Type == "Positive") %>%
  group_by(From, To) %>%
  summarize(weight = max(weight)) 

network_combined <- data_combined %>%
  select(-weight) %>%
  graph_from_data_frame(d = ., directed = F,
                        vertices = slice(Max_Count, match(unique(c(pull(., From), 
                                                                   pull(., To))),
                                                          OTU_ID)) %>%
                          mutate_if(is.factor, as.character) %>%
                          select_if(function(x) is.character(x) | is.numeric(x))) %>%
  set_vertex_attr(graph = ., name = "label", value = NA)

deg <- degree(network_combined, mode = "all")
V_size <- ifelse((log(deg)) < 3 , 2.5, (log(deg)*1.8))

layout_network <- layout_nicely(network_combined)

cluster <- tibble(Cluster = cluster_edge_betweenness(network_combined)$membership,
                  OTU_ID = V(network_combined)$name) %>%
  mutate(n = table(Cluster)[Cluster]) %>%
  mutate(Cluster = ifelse(n < 10, 99, Cluster)) %>%
  mutate(Cluster = c(seq(1, length(unique(Cluster))-1), "Others")[factor(Cluster, levels = unique(Cluster))]) %>%
  mutate(Colour = left_join(., read_csv("data/Cluster_Colour.csv"), by = "Cluster")$Colour) %>%
  mutate(Degree = deg)

write_csv(cluster, "output/SparCC/SparCC_Cluster.csv")

#### Cluster validation ####
library(vegan)
library(BiocParallel)

count_y <- datalist_Combined$Count_Data %>%
  filter(OTU_ID %in% cluster$OTU_ID) %>%
  dplyr::slice(match(cluster$OTU_ID, OTU_ID)) %>%
  select_if(is.numeric)

adonis2(count_y ~ Cluster, data = cluster)
### --> Probability that distances between groups defined by modules are obtained by chance: p < 0.001 
### (i.e. random assignment of ASVs to modules)

bplapply(seq_len(Nperm), function(x) {
  
  count_y <- datalist_Combined$Count_Data %>%
    filter(OTU_ID %in% cluster$OTU_ID) %>%
    dplyr::slice(match(cluster$OTU_ID, OTU_ID)) %>%
    select_if(is.numeric) %>%
    select_if(colSums(.) >= 8000) %>%
    t() %>%
    rrarefy(., sample = 8000) %>%
    t()
  
  perm <- adonis2(count_y ~ Cluster, data = cluster, by = "onedf")
  
})

modularity1 = cluster_edge_betweenness(network_combined)

obs <- modularity(modularity1)

Nperm = 1000

randomized.modularity <- bplapply(seq_len(Nperm), function(x){
  randomnet <- rewire(network_combined, with=each_edge(0.5))
  return(cluster_edge_betweenness(randomnet) %>% modularity())
})

lapply(randomized.modularity, function(x) x > obs) %>%
  unlist() %>% sum()/Nperm

### --> Probability that a network with comparable modularity is formed: p < 0.001
### (i.e. the modularity of the graph partitioning)
### (using random rewiring -> 50% chance for an edge to rewire terminal node)



#### Network visualization ####

png("figs/Network_Environment.png",
    res = 350, height = 250, width = 350, units = "mm", pointsize = 10)

par(mfrow = c(2,3))

plot(network_combined, vertex.size = V_size, rescale = T, layout = layout_network*0.03,
     vertex.color = get_colors_cont(vertex_attr(network_combined, "Pot_Temperature"), "RdBu", reverse = T),
     main = "Temperature")

plot(network_combined, vertex.size = V_size, rescale = T, layout = layout_network*0.03,
     vertex.color = get_colors_cont(abs(vertex_attr(network_combined, "Latitude")), "YlOrBr"),
     main = "Abs. Latitude")

plot(network_combined, vertex.size = V_size, rescale = T, layout = layout_network*0.03,
     vertex.color = get_colors_cont(vertex_attr(network_combined, "Depth"), "PuOr"),
     main = "Depth")

plot(network_combined, vertex.size = V_size, rescale = T, layout = layout_network*0.03,
     vertex.color = get_colors_cont(vertex_attr(network_combined, "Salinity"), "RdYlBu", reverse = T),
     main = "Salinity")

plot(network_combined, vertex.size = V_size, rescale = T, layout = layout_network*0.03,
     vertex.color = get_colors_cont(log(vertex_attr(network_combined, "Si")), "Reds"),
     main = "Log. Silicate")

plot(network_combined, vertex.size = V_size, rescale = T, layout = layout_network*0.03,
     vertex.color = get_colors_cont(log(vertex_attr(network_combined, "Fluorescence")), "Greens", reverse = F),
     main = "Log. Fluorescence")

dev.off()

png("figs/Network_Cluster.png",
    res = 350, height = 250, width = 350, units = "mm", pointsize = 10)

par(mfrow = c(1,1))

plot(network_combined, vertex.size = V_size, rescale = T, layout = layout_network*0.03,
     vertex.color = read_csv("output/SparCC/SparCC_Cluster.csv")$Colour,
     main = "Cluster")

dev.off()