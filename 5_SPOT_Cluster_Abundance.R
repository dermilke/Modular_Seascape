library(tidyverse)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

datalist <- import_data("../SPOT/data/Complete/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Date = lubridate::ymd(paste(Year, Month, Day, sep = "-")))

datalist_cluster <- datalist

datalist_cluster$Count_Data <- datalist %>%
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

datatable_cluster <- create_datatable(datalist_cluster, grpBy = Cluster, 
                                      otherThreshold = 0, addOthers = F, addColorScheme = F) %>%
  group_by(Date, Depth_Grp, Size_Fraction, Group) %>%
  summarize_if(is.numeric, mean) %>%
  ungroup() %>%
  group_by(Group) %>%
  mutate(Abundance_scaled = Abundance/sum(Abundance))

colours <- read_csv("output/SparCC/SparCC_Cluster.csv") %>%
  mutate(Cluster = ordered(Cluster, levels = c(as.character(seq(1, length(unique(Cluster))-1)), "Others"))) %>%
  arrange(Cluster) %>%
  select(Cluster, Colour) %>%
  distinct() %>%
  .$Colour 

datatable_cluster %>%
  filter(Depth_Grp %in% c("Surface", "DCM", "150m")) %>% 
  group_by(Size_Fraction, Group, Month, Depth_Grp) %>%
  summarize(Abundance = mean(Abundance)) %>%
  mutate(Group = ordered(Group, levels = c(as.character(1:8), "Others"), labels = c(1:8, 10))) %>%
  mutate(Depth_Grp = ordered(Depth_Grp, levels = c("Surface", "DCM", "150m", "500m", "890m"))) %>%
  filter(!is.na(Depth_Grp)) %>%
  filter(Group != "Others") %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c("D", "AE"),
                                 labels = c("0.2-1 µm", "1-80 µm"))) %>%
  ggplot(., aes(x = Month, y = Abundance * 100, fill = Group)) +
    geom_area(col = "black", size = .2) +
    scale_fill_manual(values = colours[-9]) +
    scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
    facet_grid(Depth_Grp ~ Size_Fraction) +
    labs(fill = "Cluster", y = "Abundance (%)", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1))

ggsave("figs/Network_Cluster_Abundance_Month.png", width = 7, height = 5, dpi = 300)

el_nino_table <- datatable_cluster %>%
  ungroup() %>%
  mutate(Group = ordered(Group, levels = c(as.character(1:8), "Others"), labels = c(1:8, 10))) %>%
  mutate(Depth_Grp = ordered(Depth_Grp, levels = c("Surface", "DCM", "150m", "500m", "890m"))) %>%
  filter(!is.na(Depth_Grp)) %>%
  filter(Group != "Others") %>%
  left_join(., read_csv("../SPOT/data/ENSO.csv"), by = c("Year", "Month")) %>%
  mutate(Event_Strength = paste0(Event, " - ", Strength)) %>%
  mutate(Event_Strength = ordered(Event_Strength, levels = c("El_Nino - Very Strong", "El_Nino - Moderate",
                                             "El_Nino - Weak", "NA - NA", "La_Nina - Weak", 
                                             "La_Nina - Moderate", "La_Nina - Strong"),
                                  labels = c("El Niño - Very Strong", "El Niño - Moderate",
                                             "El Niño - Weak", "No anomaly",
                                             "La Niña - Weak", "La Niña - Moderate",
                                             "La Niña - Strong"))) %>%
  
  filter(Depth_Grp == "Surface", Size_Fraction == "D") 

colors <- c('#b2182b','#f4a582','#fddbc7','white','#92c5de','#4393c3','#2166ac')

label <- lubridate::year(unique(el_nino_table$Date)) %>%
  as.character()
  
for (i in 1:(length(label)-1)) {
  tmp <- label[i]
  if (label[i+1] == tmp) label[i] <- ""
}

vline_pos <- which(label != "")+0.5

el_nino_table %>% 
  mutate(Size_Fraction = "0.22-1 µm") %>%
  ggplot(aes(x = as.factor(Date), y = Abundance*100)) +
    geom_bar(aes(fill = Group), position = "stack", stat = "identity", col = "black", size = 0.1) +
    geom_vline(xintercept = vline_pos, linetype = "dotted") +
    scale_fill_manual(values = colours[-9]) +
    scale_y_continuous(limits = c(0, 70)) +
    scale_x_discrete(breaks = as.factor(unique(el_nino_table$Date)), labels = c(label[-c(1:5)], label[1:5])) +
    facet_grid(~Size_Fraction, scale = "free_x", space = "free") +
    labs(y = "Abundance (%)", x = "") +
    theme_bw() +
      theme(axis.text.x = element_text(angle = 0, size = 9),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 

ggsave("figs/Barplot_D_Surface_Cluster.png", width = 8, height = 4, dpi = 300)
    
el_nino_table %>%
  select(Date, Event_Strength, POI) %>%
  distinct() %>%
ggplot(., aes(x = as.factor(Date), y = POI, fill = Event_Strength)) +
  geom_hline(yintercept = seq(0.5, 2, 0.5), col = c('#fddbc7', '#f4a582', '#d6604d', '#b2182b'), 
             size = .8, lty = 1) +
  geom_hline(yintercept = seq(-1.5, -0.5, 0.5), col = c('#2166ac', '#4393c3', '#92c5de'), 
             size = .8, lty = 1) +
  geom_bar(stat = "identity", col = "black", size = .2) +
  scale_fill_manual(values = colors) +
  geom_line() +
  labs(y = "Pacific Oscillation Index", x = "", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
ggsave("figs/POI_Filtered_D.png", width = 8, height = 2, dpi = 300)
