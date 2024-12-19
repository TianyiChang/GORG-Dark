library(tidyverse)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")

metag_metadata <- read.csv("metadata/metag_metat_sag_v4.csv") %>% 
    as_tibble() %>% 
    mutate(run_accessions = str_remove(run_accessions, "\\s"))

sag_metadata <- read_tsv("../sag_metadata/v3_SAG_summary_20240320.tsv")

sag_plate_region <- metag_metadata %>% 
    select(plate, ocean_province)

sag_metadata <- sag_metadata %>% 
    # add 'ocean_province'
    mutate(plate = str_replace(SAG, "^([:alpha:]+-\\d+)-.*$", "\\1")) %>% 
    left_join(sag_plate_region, by = "plate") %>% 
    # remove contaminants: AM-259-L06 and ​AM-685-K20 
    filter(SAG != "AM-259-L06") %>% 
    filter(SAG != "AM-685-K20")

write_csv(sag_metadata, "../sag_metadata/v4_SAG_summary_20240320.tsv")