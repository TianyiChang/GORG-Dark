#!/home/tchang/miniconda3/envs/shortread/bin/Rscript --vanilla

# this script predicts depth layer of each SAG based on the metag layer info and sag_norm_abund in each metag using games-howell

library(tidyverse)
library(rstatix)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")

sag_norm_abund <- read_csv("sag_abund_metag/norm_abund_sunlit_dark_sra_collab_gd_gt.csv")
metag_metadata <- read.csv("metadata/metag_metat_sag_v4.csv") %>% 
    as_tibble() %>% 
    mutate(run_accessions = str_remove(run_accessions, "\\s"))

metag_depth <- metag_metadata %>% 
    filter(!is.na(run_accessions)) %>% 
    distinct(run_accessions, depth_group)

metag_region <- metag_metadata %>% 
    filter(!is.na(run_accessions)) %>% 
    distinct(run_accessions, ocean_province)

sag_norm_abund_region <- sag_norm_abund %>% 
    left_join(metag_region, by = "run_accessions")

# handle inconsistent run name
join_norm_abund_and_metadata_var <- function(metadata_var){
    sag_norm_abund %>% 
        mutate(run_accessions = str_remove(
            run_accessions, "_S\\d{2}_HB27filtered")
            ) %>% # remove suffix for HB27filtered metag to be consistent with metag_metadata
        left_join({metadata_var}, by = "run_accessions")
}

sag_norm_abund_depth <- join_norm_abund_and_metadata_var(metag_depth)
sag_norm_abund_region <- join_norm_abund_and_metadata_var(metag_region)

#==================================================#
##  sanity check for metag lack depth layer info ##
#==================================================#

run_lack_depth_in_sag_norm_abund_depth <- filter(sag_norm_abund_depth, is.na(depth_group)) %>% 
    distinct(run_accessions)

run_lack_region_in_sag_norm_abund_region <- filter(sag_norm_abund_region, is.na(ocean_province)) %>% 
    distinct(run_accessions) # pass: no region info is missing

run_lack_depth_in_metadata <- filter(metag_metadata, run_accessions %in% run_lack_depth_in_sag_norm_abund_depth$run_accessions)

#! noted of specical cases: metag from 'Black Sea', "Ross Ice Shelf", and "Baltic Sea" lack depth layer info due to their unusal depth of aphotic zone
distinct(run_lack_depth_in_metadata, ocean_province)
# pass: All metag from 'Black Sea', "Ross Ice Shelf", and "Baltic Sea" lack depth layer info
filter(metag_metadata, ocean_province %in% c("Black Sea", "Ross Ice Shelf", "Baltic Sea")) %>% 
    filter(!is.na(depth_group))

# pass: No metag in norm_abund df is missing from the metag_metadata
metag_missing_from_metag_metadata <- filter(run_lack_depth_in_sag_norm_abund_depth, !run_accessions %in% metag_metadata$run_accessions)

#===========================================================#
##  Games-Howell test to predict depth layer of indiv SAGs ##
#===========================================================#

#! noted: there is lack of metag from the hadal zone (only 4 metag)

sag_norm_abund_no_na_depth <- sag_norm_abund_depth %>% 
    filter(!is.na(depth_group))

sag_norm_abund_region <- sag_norm_abund_region %>% 
    filter(!is.na(ocean_province)) %>% 
    filter(ocean_province != "undefined")

#?: Add a minimal mapped_count_perMread_perMbp cutoff e.g., the min mapped_count_perMread_perMbp for the SAGs with significant pair

run_games_howell <- function(df, yvar, condition, stats_out, summary_out) {

    sags_w_no_reads_mapped <- {{df}} %>%
        group_by(sag) %>% #! HERE SHOULD group_by SAG + CONDITION (e.g., depth or region)
        summarise(n=sum(mapped_count_perMread_perMbp)) %>%
        ungroup() %>%
        filter(n == "0")

    df_processed <- {{df}} %>%
        anti_join(sags_w_no_reads_mapped, by = "sag") %>% # for the test to perform properly, remove zero y-values
        filter(!is.na(mapped_count_perMread_perMbp)) %>% # three sags (AH-888-F19, AH-988-D20, AM-262-O02) lack assembly_length info, so have NAs in 'mapped_count_perMread_perMbp'
        mutate(
            sag = as.factor(sag),
            {{yvar}} := as.factor({{yvar}}),
            mapped_n_dummy =
                mapped_count_perMread_perMbp +
                rnorm(length(mapped_count_perMread_perMbp), 0.005, 0.00001) # for the test to perform properly, adding a varied small value to each row
            )
    
    #! noted: estimate = group2 - group1
    #! noted: using pvalue cutoff 0.05 to be consistent with 'frag_recruit_particle_vs_fl.r'
    if({{condition}} == "depth") {
        games_howell_summary <- df_processed %>%
            group_by(sag) %>%
            games_howell_test(mapped_n_dummy ~ depth_group, conf.level = 0.95, detailed = TRUE)
    } else if ({{condition}} == "region") {
        games_howell_summary <- df_processed %>%
            group_by(sag) %>%
            games_howell_test(mapped_n_dummy ~ ocean_province, conf.level = 0.95, detailed = TRUE)
    } else {
       print("Condition is not existed")
    }

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

    # To find the highest hierarchy of depth/region, the parent should not be found in child
    major_group_assigned <- games_howell_signf_child %>% 
        left_join(games_howell_signf_parent, by = "sag") %>% 
        group_by(sag) %>% 
        filter(!parent %in% child) %>% 
        distinct(sag, parent) %>% 
        mutate(multiple_major_group = ifelse(
                n() > 1, "yes", "no"))        
    
    # sort the groups
    if({{condition}} == "depth"){
        major_group_sorted <- major_group_assigned %>% 
            mutate(
                parent = factor(parent, levels = c(
                    "epi", "meso", "bathy", "abysso", "hadal"))
            ) %>% 
            arrange(parent)
    } else if ({{condition}} == "region") {
        major_group_sorted <- major_group_assigned %>% 
            mutate(
                parent = factor(parent, levels = c(
                    "North Pacific Ocean", "South Pacific Ocean",
                    "North Atlantic Ocean", "South Atlantic Ocean",
                    "Indian Ocean", "Arctic Ocean",
                    "Southern Ocean", "Ross Ice Shelf",
                    "Mediterranean Sean", "Red Sea",
                    "Baltic Sea", "Black Sea"
                    )) # sort the depth layer
            ) %>% 
            arrange(parent)
    } else {
       print("Condition is not existed")
    }

    major_group_4_indiv_sags <- major_group_sorted %>% 
        # combine multiple main groups
        mutate(major_group_index = str_c(
            "major_group_", c(1:n())
        )) %>% 
        ungroup() %>% 
        pivot_wider(
            names_from = major_group_index,
            values_from = parent
            ) %>% 
        unite(major_group, starts_with("major_group_"), remove = TRUE, sep = " - ") %>% 
        select(sag, major_group, multiple_major_group) %>% 
        mutate(major_group = str_remove(major_group, " - NA -.*| - NA$"))

    # export final results
    write.csv(games_howell_summary, str_c("sag_depth_region_preference/", {{stats_out}}))
    write_csv(major_group_4_indiv_sags, str_c("sag_depth_region_preference/", {{summary_out}}))

    # return(
    #     df_list = list(
    #         games_howell_summary = games_howell_summary,
    #         major_depth_4_indiv_sags = major_depth_4_indiv_sags
    #     )
    # )

}

print("Running games-howell test for depth preference")
run_games_howell(
    sag_norm_abund_no_na_depth, depth_group,
    "depth", "sunlit_dark_sag_depth_games_howell_summary.csv",
    "sunlit_dark_major_depth_4_indiv_sags.csv"
    )

print("Running games-howell test for region preference")
run_games_howell(
    sag_norm_abund_region, ocean_province,
    "region", "sunlit_dark_sag_region_games_howell_summary.csv",
    "sunlit_dark_major_region_4_indiv_sags.csv"
    )
