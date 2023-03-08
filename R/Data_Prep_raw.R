library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")

datalist_Atlantic <- import_data("~/PhD/Data_Storage/ASV_Tabs/Atlantic/V4V5_Primerset/Complete/", kingdom = "Chloroplast", 
                                 abundance_filter = F, min_counts = 0) %>%
  mutate_meta_datalist(Station = as.character(Station)) %>%
  mutate_meta_datalist(Ocean = "Atlantic")

datalist_Pacific <- import_data("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", kingdom = "Chloroplast", 
                                 abundance_filter = F, min_counts = 0)%>%
  mutate_meta_datalist(Ocean = "Pacific")

datalist_Combined <- combine_data(
  datalist_Pacific %>% 
    filter_station_datalist(Latitude > max(c(min(datalist_Atlantic$Meta_Data$Latitude), 
                                             min(datalist_Pacific$Meta_Data$Latitude))) &
                            Latitude < min(c(max(datalist_Atlantic$Meta_Data$Latitude), 
                                             max(datalist_Pacific$Meta_Data$Latitude)))) %>%
    filter_station_datalist(Depth <= 200),

  datalist_Atlantic %>% 
    filter_station_datalist(Latitude > max(c(min(datalist_Atlantic$Meta_Data$Latitude), 
                                             min(datalist_Pacific$Meta_Data$Latitude))) &
                            Latitude < min(c(max(datalist_Atlantic$Meta_Data$Latitude), 
                                             max(datalist_Pacific$Meta_Data$Latitude))))
)

query <- tibble(OTU_ID = datalist_Combined$Count_Data$OTU_ID, 
                Sequence = map_chr(datalist_Combined$Count_Data$OTU_ID, 
                                   function(x) system(paste0("grep -A 1 ", x, " ", 
                                                      "~/PhD/Data_Storage/Fasta/Complete/V4V5_Primerset/Prok/Sequences_Complete.fasta"), 
                                                      intern = T)[2])) %>%
  mutate(OTU_ID = paste0(">", OTU_ID))

is_na_query <- query %>%
  filter(is.na(Sequence)) %>%
  mutate(OTU_ID = str_replace(OTU_ID, pattern = "^>", replacement = ""))

from_atlantic_fasta <- is_na_query %>%
  mutate(Sequence = 
    map_chr(is_na_query$OTU_ID, 
          function(x) system(paste0("grep -A 1 ", x, " ", 
                                    "~/PhD/Data_Storage/Fasta/Atlantic/V4V5_Primerset/Prok/Sequences_Atlantic_Complete.fasta"), 
                             intern = T)[2])
  )

also_Pacific <- from_atlantic_fasta %>%
  filter(is.na(Sequence)) %>%
  mutate(Sequence = 
           map_chr(OTU_ID, 
                   function(x) system(paste0("grep -A 1 ", x, " ", 
                                             "~/PhD/Data_Storage/Fasta/Pacific/V4V5_Primerset/Prok/Sequences_Pacific_Complete.fasta"), 
                                      intern = T)[2])
  )

second_Pacific <- also_Pacific %>%
  filter(is.na(Sequence)) %>%
  mutate(Sequence = 
           map_chr(OTU_ID, 
                   function(x) system(paste0("grep -A 1 ", x, " ", 
                                             "~/PhD/Data_Storage/Fasta/Pacific/V4V5_Primerset/Euk/Pacific_Pool_6.fasta"), 
                                      intern = T)[2])
  )

final <- query %>%
  filter(!is.na(Sequence)) %>%
  bind_rows(., filter(from_atlantic_fasta, !is.na(Sequence)), also_Pacific) 

do.call(rbind, lapply(seq(nrow(final)), function(i) t(final[i, ]))) %>%
  write.table(., file = "~/PhD/Projects/Comparison_Paper/data/Fasta/Chloroplast/Chloroplast_Sequences.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)

