#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

library(tidyverse)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/related_studies/omd")

mags <- read_tsv("genomes-usage.tsv") %>%
    mutate(
        dataset=str_replace(genome, "_.*", ""),
        file=str_replace(genome, "$", ".fa")
        ) %>%
    filter(dataset != "GORG" & dataset != "MARD") %>% 
    filter(bacteria == "TRUE" | archaea == "TRUE") 

write_lines(mags$file, "mag_list.txt")