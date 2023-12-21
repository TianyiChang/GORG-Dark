#!/home/tchang/miniconda3/envs/shortread/bin/Rscript --vanilla

# this script predicts depth layer of each SAG based on the metag layer info and sag_norm_abund in each metag using games-howell

library(tidyverse)
library(rstatix)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")

sag_norm_abund <- read_csv("sag_abund_metag/norm_abund_dark_collab_metag_gdv3.csv")
metag_metadata <- read_csv("metadata/metag_metat_sag_v3.csv")

metag_depth <- metag_metadata %>% 
    filter(!is.na(run_accessions)) %>% 
    distinct(run_accessions, depth_group)

sag_norm_abund_depth <- sag_norm_abund %>% 
    mutate(run_accessions = str_remove(
        run_accessions, "_S\\d{2}_HB27filtered")
        ) %>% # remove suffix for HB27filtered metag to be consistent with metag_metadata
    left_join(metag_depth, by = "run_accessions")

#==================================================#
##  sanity check for metag lack depth layer info ##
#==================================================#

run_lack_depth_in_sag_norm_abund_depth <- filter(sag_norm_abund_depth, is.na(depth_group)) %>% 
    distinct(run_accessions)

run_lack_depth_in_metadata <- filter(metag_metadata, run_accessions %in% run_lack_depth_in_sag_norm_abund_depth$run_accessions)

#! noted of specical cases: metag from 'Black Sea', "Ross Ice Shelf", and "Baltic Sea" lack depth layer info due to their unusal depth of aphotic zone
distinct(run_lack_depth_in_metadata, ocean_province)
# pass: All metag from 'Black Sea', "Ross Ice Shelf", and "Baltic Sea" lack depth layer info
filter(metag_metadata, ocean_province %in% c("Black Sea", "Ross Ice Shelf", "Baltic Sea"))

# pass: No metag in norm_abund df is missing from the metag_metadata
metag_missing_from_metag_metadata <- filter(run_lack_depth_in_sag_norm_abund_depth, !run_accessions %in% metag_metadata$run_accessions)

#===========================================================#
##  Games-Howell test to predict depth layer of indiv SAGs ##
#===========================================================#

#! noted: the df contains epi metag from the Indian Ocean, as for now, there is no norm abund from epi zone of other oceans
#! noted: there is lack of metag from the hadal zone (only 2 metag were collected, and were discarded in frag_recruit_sag_norm_abund.r L32)

sag_norm_abund_depth_no_epi <- sag_norm_abund_depth %>% 
    filter(depth_group != "epi")

sag_norm_abund_no_na_depth <- sag_norm_abund_depth_no_epi %>% 
    filter(!is.na(depth_group))

sags_w_no_reads_mapped <- sag_norm_abund_no_na_depth %>%
    group_by(sag) %>%
    summarise(n=sum(mapped_count_perMread_perMbp)) %>%
    ungroup() %>%
    filter(n == "0")

sag_norm_abund_no_na_depth_no_zeros <-
    sag_norm_abund_no_na_depth %>%
    anti_join(sags_w_no_reads_mapped, by = "sag") %>% # for the test to perform properly, remove zero y-values
    mutate(
        sag = as.factor(sag),
        depth_group = as.factor(depth_group),
        mapped_n_dummy =
            mapped_count_perMread_perMbp +
            rnorm(length(mapped_count_perMread_perMbp), 0.005, 0.00001) # for the test to perform properly, adding a varied small value to each row
        )

#! noted: estimate = group2 - group1
#! noted: using pvalue cutoff 0.05 to be consistent with 'frag_recruit_particle_vs_fl.r'
games_howell_summary <- sag_norm_abund_no_na_depth_no_zeros %>%
    group_by(sag) %>%
    games_howell_test(mapped_n_dummy ~ depth_group, conf.level = 0.95, detailed = TRUE) %>%
    left_join(sag_norm_abund_no_na_depth_no_zeros, by = "sag")

games_howell_signf <- games_howell_summary %>% 
    filter(p.adj < 0.05) %>% 
    # specify child and parent based on estimate
    group_by(sag) %>% 
    mutate(
        child = ifelse(estimate > 0, group1, group2),
        parent = ifelse(estimate > 0, group2, group1)
    ) %>% 
    ungroup()

games_howell_signf_parent <- games_howell_signf %>% 
    distinct(sag, parent)

games_howell_signf_child <- games_howell_signf %>% 
    distinct(sag, child)

# To find the highest hierarchy of depth_group, the parent should not be found in child
major_depth_4_indiv_sags <- games_howell_signf_child %>% 
    left_join(games_howell_signf_parent, by = "sag") %>% 
    group_by(sag) %>% 
    filter(!parent %in% child) %>% 
    distinct(sag, parent) %>% 
    mutate(multiple_major_depth = ifelse(
        n() > 1, "yes", "no"
    ))

# export final results
write.csv(games_howell_summary, "sag_depth_preference/sag_depth_games_howell_summary.csv")
write_csv(major_depth_4_indiv_sags, "sag_depth_preference/major_depth_4_indiv_sags.csv")