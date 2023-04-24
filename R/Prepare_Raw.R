read_taxonomy <- function(file, kingdom, DB, fillEmpty = T) {
  
  fill_empty <- function(x, entity) {
    
    for (i in (3:ncol(select_if(x, is.character)))) {
      
      x[is.na(x[,i]),i] <- x[is.na(x[,i]),i-1] %>%
        pull(.) %>%
        paste(entity, ., sep = " ")
      
      x[,i] <- gsub(paste(entity, " ", entity, " ", sep = ""), paste(entity, " ", sep = ""), pull(x[,i]))
      
    }
    
    return(x)
    
  }
  
  if (kingdom == "Prok") {
    
    taxonomy_prok <- suppressWarnings(
      suppressMessages(read_tsv(file)) %>%
        dplyr::slice(c(grep("d__Bacteria",  #"d__Bacteria", 
                            .$taxonomy, value = F),
                       grep("d__Archaea", .$taxonomy, value = F))) %>%
        separate(., 2, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
        mutate_if(is.character, 
                  str_replace_all, pattern = " *[a-z]__", #"D_[0-9]*__"
                  replacement = "") %>%
        dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
    )
    
    taxonomy_chloroplast <- suppressWarnings(
      suppressMessages(read_tsv(file)) %>%
        dplyr::slice(grep("Kingdom.Eukaryota", .$taxonomy, value = F)) %>%
        separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Subclass","Order","Suborder", "Family", "Genus", "Species"), sep = ";") %>%
        mutate_if(is.character, 
                  str_replace_all, pattern = ".*\\.", replacement = "") %>%
        dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
        with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
    )
    
    taxonomy <- list(Prok = taxonomy_prok, Chloroplast = taxonomy_chloroplast)
    
  } else if (kingdom == "Euk") {
    
    if (DB == "SILVA132") {
      
      taxonomy <- suppressWarnings(
        suppressMessages(read_tsv(file)) %>%
          separate(., 2, c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = ";") %>%
          mutate_if(is.character, 
                    str_replace_all, pattern = "D_[0-9]*__", replacement = "") %>%
          dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
          with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
      )
      
    } else if (DB == "PR2") {
      
      taxonomy <- suppressWarnings(
        suppressMessages(read_tsv(file)) %>%
          separate(., 2, c("Kingdom","Supergroup","Phylum","Class","Order", "Family", "Genus", "Species"), sep = ";") %>%
          mutate_if(is.character, 
                    str_replace_all, pattern = "[a-z]*_", replacement = "") %>%
          dplyr::rename(., 'OTU_ID' = '#OTUID') %>%
          with(., if (fillEmpty) {fill_empty(., "unknown")} else {.})
      )
      
    }
  }
  
  return(taxonomy)
  
}

prepare_raw <- function(file_ASV, file_Meta, confidence_lvl = 0.8, kingdom = "Prok", DB = "SILVA132", fill = T) {
  
  if (kingdom == "Euk") {
    
    taxonomy <- read_taxonomy(paste(file_ASV, "Raw/", kingdom, "/taxonomy-", DB, ".tsv", sep = ""), kingdom, DB, fill)
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-18S-seqs.with-", DB, "-tax.tsv", sep = ""),
                                              delim = "\t", skip = 1)) %>%
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
      right_join(taxonomy, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data, paste(file_ASV, "Processed/Euk/Full_Euk_Count.tsv", sep = ""))
    
    Meta_Data <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      dplyr::slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data, paste(file_ASV, "Meta_Data/Euk/Meta_Data.tsv", sep = ""))
    
  } else {
    
    taxonomy <- read_taxonomy(paste(file_ASV, "Raw/Prok/taxonomy.tsv", sep = ""), kingdom, fillEmpty = fill)
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-16S-seqs.with-tax.tsv", sep = ""), 
                                              delim = "\t", skip = 1)) %>%
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
      left_join(taxonomy$Prok, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data, paste(file_ASV, "Processed/Prok/Full_Prok_Count.tsv", sep = ""))
    
    Meta_Data <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      dplyr::slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data, paste(file_ASV, "Meta_Data/Prok/Meta_Data.tsv", sep = ""))
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Raw/", kingdom, "/all-16S-seqs.with-tax.tsv", sep = ""), 
                                              delim = "\t", skip = 1)) %>%
      dplyr::rename(., 'OTU_ID' = '#OTU ID') %>%
      left_join(taxonomy$Chloroplast, ., by = 'OTU_ID') %>%
      filter(confidence >= confidence_lvl) %>%                                                           
      select(-confidence) %>%
      select_if(!colnames(.) == "taxonomy") %>%
      filter(rowSums(select_if(.,is.numeric)) != 0)
    
    write_tsv(Count_Data, paste(file_ASV, "Processed/Chloroplast/Full_Chloroplast_Count.tsv", sep = ""))
    
    Meta_Data <- suppressMessages(read_delim(file_Meta, del = "\t")) %>% 
      dplyr::slice(match(names(select_if(Count_Data, is.numeric)), .$Sample_ID))
    
    write_tsv(Meta_Data, paste(file_ASV, "Meta_Data/Chloroplast/Meta_Data.tsv", sep = ""))
    
  }
  
  cat("Raw tables converted into ", paste(file_ASV, "\n", sep = ""))
  
}

data_select <- function(file_ASV, kingdom = "Prok") {
  
  Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Processed/", kingdom, "/Full_", kingdom,"_Count.tsv", sep = ""), 
                                            del = "\t")) %>%
    select(-grep("[Mm]ock", names(.))) %>%
    select(-grep("[Bb]lank", names(.))) %>%
    select(-grep("staggered-16s", names(.))) %>%
    select(-grep("even-16s", names(.))) %>%
    select(-grep("NC", names(.))) 
  
  Meta_Data <- suppressMessages(read_delim(paste(file_ASV, "Meta_Data/", kingdom, "/Meta_Data.tsv", sep = ""), 
                                           del = "\t"))
  
  return(list(Meta_Data = Meta_Data, 
              Count_Data = Count_Data))
  
}