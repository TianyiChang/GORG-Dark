#!/home/tchang/miniconda3/envs/shortread/bin/Rscript --vanilla

library(tidyverse)
library(ape)
library(vegan)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")

#===========================#
## !Remove epi metag data ##
#===========================#

sag_norm_abund_ori <- read_csv("sag_abund_metag/norm_abund_dark_collab_metag_gdv3.csv")
metag_metadata <- read_csv("metadata/metag_metat_sag_v3.csv")

metag_depth <- metag_metadata %>% 
    filter(!is.na(run_accessions)) %>% 
    distinct(run_accessions, depth_group)

# handle inconsistent run name
sag_norm_abund <- sag_norm_abund_ori %>% 
    mutate(run_accessions = str_remove(
        run_accessions, "_S\\d{2}_HB27filtered")
        ) %>% # remove suffix for HB27filtered metag to be consistent with metag_metadata
    left_join(metag_depth, by = "run_accessions") %>% 
    filter(depth_group != "epi") %>% 
    select(sag, run_accessions, mapped_count_perMread_perMbp)

#============#
## Run PCoA ##
#============#

# multivariate analyses
# points: individ. sags
# df: columns (samples), rows (sags)

sag_norm_abund_zero <- sag_norm_abund %>%
    group_by(sag) %>%
    summarise(sum = sum(mapped_count_perMread_perMbp)) %>%
    filter(sum == "0")

non_zero_df <- sag_norm_abund %>%
    anti_join(sag_norm_abund_zero, by = "sag") %>%
    pivot_wider(names_from = run_accessions,
        values_from = mapped_count_perMread_perMbp) %>%
    column_to_rownames(var = "sag")

#PCoA from ape
bray_pcoa <- non_zero_df %>%
    vegdist(method = "bray") %>%
    pcoa()

# get eigenvalues for the two major axeses
relative_eig <- tibble(
    relative_eig1 = round(bray_pcoa$values$Relative_eig[1], 2),
    relative_eig2 = round(bray_pcoa$values$Relative_eig[2], 2)
)

# get PC coord
pc_1_2 <- bray_pcoa$vectors[,c(1,2)] %>%
    as_tibble()

# combine to get final df
df <- sag_norm_abund %>%
    anti_join(sag_norm_abund_zero, by = "sag") %>%
    distinct(sag) %>%
    mutate(plate = str_replace(
            sag, "(.*)-.*", "\\1")) 

df_2 <- bind_cols(df, pc_1_2) %>% 
    rename(PC_1 = Axis.1, PC_2 = Axis.2)

write_csv(df_2, "pcoa/sag_pc1_2_values.csv")