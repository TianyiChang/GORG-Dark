#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

library(tidyverse)
library(furrr)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/related_studies/omd/omd_m_dark")

plan(multisession, workers = 30)

get_depth <- function(path2file) {

    df <- read_csv({{path2file}})
    
    depth_col <- df %>%
        select(matches("(?i)depth")) %>%
        mutate_if((is.character), parse_number)

    out <- tibble(accession = df$accession, depth_col)

    return(out)
}

combined <- list.files(
    path = "metad/df",
    full.names = TRUE,
    pattern = ".csv") %>% 
    future_map_dfr(~get_depth(.))

# note: only two selected columns ('depth' and 'Depth') for the current data
acc_depth <- combined %>% 
    mutate(depth = ifelse(
        is.na(Depth), depth, Depth)) %>% 
    select(-Depth)


omd_m <- tibble(file = list.files(
    path = "omd_mags", pattern = ".fa")) %>% 
    mutate(accession = str_replace(
        file, "^[^_]+_([^_]+)_.*$", "\\1")) %>% 
    left_join(acc_depth, by = "accession")

# keep only mags from >200m samples
omd_m_filterd <- omd_m %>% 
    filter(depth > 200)

# export dark ocean omd file list
write_lines(omd_m_filterd$file, "omd_dark.txt")
