#### prepare Fasta for SINA alignment ####

library(tidyverse)

source("R/Import_Data.R")
source("../Jonatan_Project/R/Datalist_Wrangling_Functions.R")
source("R/Import_SparCC_Network.R")

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  filter_taxa_datalist(Family != "Mitochondria")

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  filter_taxa_datalist(Family != "Mitochondria")

prok_fasta_Atlantic <- seqinr::read.fasta("data/Atlantic/Fasta/Prok/Prok_Sequences.fasta")
prok_fasta_Pacific <- seqinr::read.fasta("data/Pacific/Fasta/Prok/Prok_Sequences.fasta")

ind <- NULL

inds <- map(datalist_Atlantic$Count_Data$OTU_ID, function(x) {
  ind <- c(ind, match(x, names(prok_fasta_Atlantic))[1])
}) %>%
  as_vector()

prok_fasta_subset_Atlantic_1 <- prok_fasta_Atlantic[inds[1:1000]]
prok_fasta_subset_Atlantic_2 <- prok_fasta_Atlantic[inds[1001:length(inds)]]

prok_fasta_subset_Pacific_1 <- prok_fasta_Pacific[datalist_Pacific$Count_Data$OTU_ID[1:1000]]
prok_fasta_subset_Pacific_2 <- prok_fasta_Pacific[datalist_Pacific$Count_Data$OTU_ID[1001:nrow(datalist_Pacific$Count_Data)]]

seqinr::write.fasta(sequences = prok_fasta_subset_Atlantic_1, names = names(prok_fasta_subset_Atlantic_1),
                    file.out = "data/Atlantic/Fasta/Prok/Prok_Atlantic_Subset_1.fasta")

seqinr::write.fasta(sequences = prok_fasta_subset_Atlantic_2, names = names(prok_fasta_subset_Atlantic_2),
                    file.out = "data/Atlantic/Fasta/Prok/Prok_Atlantic_Subset_2.fasta")

seqinr::write.fasta(sequences = prok_fasta_subset_Pacific_1, names = names(prok_fasta_subset_Pacific_1),
                    file.out = "data/Pacific/Fasta/Prok/Prok_Pacific_Subset_1.fasta")

seqinr::write.fasta(sequences = prok_fasta_subset_Pacific_2, names = names(prok_fasta_subset_Pacific_2),
                    file.out = "data/Pacific/Fasta/Prok/Prok_Pacific_Subset_2.fasta")

### Then use online-version of SINA aligner in standard settings ####

#### Create Tree with FastTree ####

system(paste0("~/PhD/Statistics/FastTree/FastTree -gtr -nt < ",
              "data/Atlantic/Aligned/Prok/SINA_Aligned_Prok_Atlantic.fasta > ",
              "data/Atlantic/Tree/Prok/Atlantic_Prok_FastTree.tree"))

system(paste0("~/PhD/Statistics/FastTree/FastTree -gtr -nt < ",
              "data/Pacific/Aligned/Prok/SINA_Aligned_Prok_Pacific.fasta > ",
              "data/Pacific/Tree/Prok/Pacific_Prok_FastTree.tree"))

