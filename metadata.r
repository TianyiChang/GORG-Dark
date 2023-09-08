#!/home/tchang/miniconda3/envs/rgdal/bin/Rscript

#todo: review and simplify the below code

#! this script aims to collect public metag samples covering geographic locations widely and evenly,
#! and the ones located within the same location cluster with metatranscriptome or surface-metag samples will have higher prority

# cluster collected samples based on longtitude and latitude
# reference: https://gis.stackexchange.com/questions/17638/clustering-spatial-data-in-r

#! note 202309: adjust the aphotic zone for Black Sea plates to >45m in "addit_seas.csv" (now is "addit_seas_v2.csv")
#! ref: https://www.frontiersin.org/articles/10.3389/fmars.2017.00090/full
#! ref: https://www.mdpi.com/2673-1924/1/4/18

library(rgdal)
library(geosphere)
library(tidyverse)
library(RColorBrewer)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/221007")

sag_meta_location <- read_csv("overview_gorg_dark_v1_20221104.csv") %>%
    mutate(
        group = "sag",
        plate = str_replace(sample_id, "^(\\w+-\\w+)-.*$", "\\1"),
        depth = as.numeric(
            str_replace(depth, "\\D+", "")
        ),
        depth_group = case_when(
            depth <= 200 ~ "epi",
            depth > 200 & depth <= 1000 ~ "meso",
            depth > 1000 & depth <= 4000 ~ "bathy",
            depth > 4000 & depth <= 6000 ~ "abysso",
            depth > 6000 ~ "hadal"
            )
        ) %>%
    distinct(plate, .keep_all = TRUE) %>%
    select(latitude, longitude, group, depth, depth_group, plate)

# inclue Black Sea plate info manually
black_sea_location <- data_frame(
    latitude = c(43.5320233, 43.5320233),
    longitude = c(32.5151083, 32.5151083),
    group = c("sag", "sag"),
    depth = c(94.7, 104.6),
    depth_group = c(NA, NA),
    plate = c("AM-685", "AM-689")
)

sag_meta_location <- bind_rows(
    sag_meta_location, black_sea_location
)

# tara_metat_sample <- read_csv("metat_tara_prok_vir.csv") %>%
#     select(sample_accessions, run_accessions) %>%
#     separate_rows(run_accessions, sep = ";") # this will create new rows

# function for metadata downloaded using sra-run-selector
get_info_from_run_selector <- function(file) {

    read_csv(file = file) %>%
        select(run_accessions = Run, longitude = Longitude_Start,
            latitude = Latitude_Start, depth = Depth,
            sample_accessions = BioSample, BioProject) %>%
        mutate(
            group = case_when(
                BioProject == "PRJNA289734" ~ "uncollected",
                BioProject == "PRJNA305355" ~ "uncollected",
                TRUE ~ "metat"
            ),
            depth = as.numeric(str_replace(depth, "^\\w+-(\\w+)$", "\\1")),
            depth_group = case_when(
                depth <= 200 ~ "epi",
                depth > 200 & depth <= 1000 ~ "meso",
                depth > 1000 & depth <= 4000 ~ "bathy",
                depth > 4000 & depth <= 6000 ~ "abysso",
                depth > 6000 ~ "hadal"
                )
            ) %>%
        select(-BioProject) %>%
        separate_rows(run_accessions, sep = ";")
    
}

# apply the function on metatransciptomes
metat_tara_prok_vir <- list.files(pattern = "run_selector_metat",
    full.names = TRUE) %>%
    map_dfr(get_info_from_run_selector)

# apply the function on Mediterranean (0.22-5) and Red sea (size frac: 0.1-1.2) metagenomes
other_metag <- list.files(pattern = "run_selector_metag",
    full.names = TRUE) %>%
    map_dfr(get_info_from_run_selector) %>%
    mutate(metat_same_sample = ifelse(
        sample_accessions %in% metat_tara_prok_vir$sample_accessions,
            1, 0))

omd_metag <- read_csv("metag_all_prok_vir_metadata_omd.csv") %>%
    select(sample_name, size_fraction, latitude, longitude, depth, ocean_province) %>%
    mutate(depth = as.numeric(
        str_replace(depth, "^\\w+-(\\w+)$", "\\1")
        )) # keep lower depth bound for cases like 10-100

get_info_metag <- function(file) {

    omd_select <- read_csv(file = file) %>%
        left_join(omd_metag, by = "sample_name") %>%
        filter(size_fraction != "<-0.22" & size_fraction != "0.1-0.22") %>% #! Noted the different size fractions
        mutate(
            group = "uncollected",
            sample_accessions = str_replace(sample_name, "^\\w+_(\\w+)_.*$", "\\1"),
            metat_same_sample = ifelse(
                sample_accessions %in% metat_tara_prok_vir$sample_accessions,
                1, 0),
            depth_group = case_when(
                depth <= 200 ~ "epi",
                depth > 200 & depth <= 1000 ~ "meso",
                depth > 1000 & depth <= 4000 ~ "bathy",
                depth > 4000 & depth <= 6000 ~ "abysso",
                depth > 6000 ~ "hadal"
            )
            ) %>%
        select(run_accessions, longitude, latitude,
            metat_same_sample, group, depth, depth_group) %>%
        separate_rows(run_accessions, sep = ";")
    
    return(omd_select)

}

metag_tara_prok_vir <- get_info_metag("metag_tara_prok_vir.csv")
metag_non_tara_prok_vir <- get_info_metag("metag_non_tara_prok_vir.csv")

#! remove tara data: 9 out of the 79 collected runs had wrong long lat data
#! change redundant run_accession 'SRR4028224' to 'SRR4028175' for metagenome 'Moc_st21_Atlantic_776'

collected_meta <- read_csv("maria_collected_221031.csv") %>%
    filter(Expedition != "TARA") %>%
    mutate(
        run_accessions = replace(run_accessions,
            metagenome == "Moc_st21_Atlantic_776",
            "SRR4028175"),
        group = "collected",
        metat_same_sample = ifelse(
            run_accessions %in% metat_tara_prok_vir$run_accessions,
            1, 0),
        depth_group = case_when(
            depth <= 200 ~ "epi",
            depth > 200 & depth <= 1000 ~ "meso",
            depth > 1000 & depth <= 4000 ~ "bathy",
            depth > 4000 & depth <= 6000 ~ "abysso",
            depth > 6000 ~ "hadal"
        )
    ) %>%
    select(run_accessions, longitude, latitude, metat_same_sample, group, depth, depth_group)

addit_collab <- read_csv("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/metadata/addit_collab.csv") %>%
    rename(run_accessions = metag, group = dataset)

addit_seas <- read_csv("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/metadata/addit_seas_v2.csv") %>%
    rename(run_accessions = run, group = dataset) %>%
    select(-depth_group)

df <- bind_rows(
    collected_meta,
    metat_tara_prok_vir,
    metag_non_tara_prok_vir,
    metag_tara_prok_vir,
    other_metag,
    sag_meta_location
) %>%
    filter(!is.na(longitude) | !is.na(latitude))

df_addit <- bind_rows(
    df, addit_collab, addit_seas
    ) %>%
    mutate(ocean_province =
        case_when(
            plate == "AM-685" | plate == "AM-689" ~ "Black Sea",
            TRUE ~ ocean_province
        )) 

df_dark <- filter(df, depth_group != "epi")

get_cluster <- function(data_frame) {

    # convert data to a SpatialPointsDataFrame object
    xy <- SpatialPointsDataFrame(
        matrix(c(data_frame$longitude, data_frame$latitude), ncol = 2), data.frame(ID = seq_along(1:nrow(data_frame))),
        proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

    # use the distm function to generate a geodesic distance matrix in meters
    mdist <- distm(xy)

    # cluster all points using a hierarchical clustering approach
    hc <- hclust(as.dist(mdist), method="complete")

    # define the distance threshold, in this case 277.8 KM = 150 Nautical Miles
    d=277800

    # define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
    data_frame$clust <- cutree(hc, h=d)

    return(data_frame)

}

cluster_all <- get_cluster(df)
cluster_dark <- get_cluster(df_dark)
cluster_df_addit <- get_cluster(df_addit)

# check the coverage for collected and uncollected
total_uniq_clust <- cluster_dark %>%
    distinct(group, clust) %>%
    mutate(count = 1) %>%
    group_by(group) %>%
    summarise(total_uniq_clust = sum(count))

# label interested avail, and select a representative each cluster, for dark and metag
# ocean provinces are defined

full_df <- cluster_all %>%
    group_by(clust) %>%
    mutate(
        metat_same_clust = ifelse("metat" %in% group,
            1, 0),
        surf_same_clust = ifelse("epi" %in% depth_group,
            1, 0),
        sag_same_clust = ifelse("sag" %in% group,
            1, 0),
        ocean_province = case_when(
            latitude >= 60 ~ "Arctic Ocean",
            latitude <= -60 ~ "Southern Ocean",
            latitude > 12 & latitude < 30 & longitude > 32 & longitude < 45 ~ "Red Sea",
            latitude > 30 & latitude < 49 & longitude > -6.5 & longitude < 44 ~ "Mediterranean Sea",
            (latitude > 0 & latitude < 60 & longitude > 105 & longitude <= 180) |
                (latitude > 0 & latitude < 60 & longitude >= -180 & longitude < -97.5) |
                    (latitude > 0 & latitude < 8 & longitude >= -97.5 & longitude < -75) ~ "North Pacific Ocean",
            (latitude < 0 & latitude > -60 & longitude > 153 & longitude <= 180) |
                (latitude < 0 & latitude > -60 & longitude >= -180 & longitude < -60)  ~ "South Pacific Ocean",
            (latitude > 0 & longitude > -75 & longitude < -6.5) |
                (latitude > 17 & latitude < 30 & longitude > -100 & longitude < -75) ~ "North Atlantic Ocean",
            latitude < 0 & longitude > -68 & longitude < 22.5 ~ "South Atlantic Ocean",
            (latitude < 0 & latitude > -60 & longitude > 22.5 & longitude < 153) |
                (latitude < 30 & latitude > 0 & longitude > 45 & longitude < 105) ~ "Indian Ocean",
            TRUE ~ "undefined"
        )
    ) %>%
    select(-sample_accessions) %>%
    ungroup()

full_df_addit <- cluster_df_addit %>%
    group_by(clust) %>%
    mutate(
        metat_same_clust = ifelse("metat" %in% group,
            1, 0),
        surf_same_clust = ifelse("epi" %in% depth_group,
            1, 0),
        sag_same_clust = ifelse("sag" %in% group,
            1, 0),
        ocean_province_2 = case_when(
            latitude >= 60 ~ "Arctic Ocean",
            latitude <= -60 ~ "Southern Ocean",
            latitude > 12 & latitude < 30 & longitude > 32 & longitude < 45 ~ "Red Sea",
            latitude > 30 & latitude < 49 & longitude > -6.5 & longitude < 44 ~ "Mediterranean Sea",
            (latitude > 0 & latitude < 60 & longitude > 105 & longitude <= 180) |
                (latitude > 0 & latitude < 60 & longitude >= -180 & longitude < -97.5) |
                    (latitude > 0 & latitude < 8 & longitude >= -97.5 & longitude < -75) ~ "North Pacific Ocean",
            (latitude < 0 & latitude > -60 & longitude > 153 & longitude <= 180) |
                (latitude < 0 & latitude > -60 & longitude >= -180 & longitude < -60)  ~ "South Pacific Ocean",
            (latitude > 0 & longitude > -75 & longitude < -6.5) |
                (latitude > 17 & latitude < 30 & longitude > -100 & longitude < -75) ~ "North Atlantic Ocean",
            latitude < 0 & longitude > -68 & longitude < 22.5 ~ "South Atlantic Ocean",
            (latitude < 0 & latitude > -60 & longitude > 22.5 & longitude < 153) |
                (latitude < 30 & latitude > 0 & longitude > 45 & longitude < 105) ~ "Indian Ocean",
            (latitude == "57.32" & longitude == "20.05") ~ "Baltic Sea",
            TRUE ~ "undefined"
        )
    ) %>%
    ungroup() %>%
    mutate(
        ocean_province = ifelse(
            is.na(ocean_province), ocean_province_2, ocean_province
        )
    ) %>%
    select(-sample_accessions, -ocean_province_2)

all_dark <- full_df %>%
    filter(depth_group != "epi" & (group == "collected" | group == "uncollected")) %>%
    arrange(-metat_same_sample, -metat_same_clust, -surf_same_clust,
        group, .by_group = TRUE)


run_represent <- all_dark %>%
    select(run_accessions, clust, metat_same_sample,
        metat_same_clust, surf_same_clust,
        longitude, latitude, depth_group) %>%
    group_by(clust) %>%
    slice_head(n = 1) %>%
    ungroup()

# save the list for downloading sra runs
setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/sra")
# write_lines(run_represent$run_accessions, "selected_sra_run.txt")
# write_lines(all_dark$run_accessions, "all_dark_sra_run.txt")

# below is newly added code 202309 for frag_recruit.smk
dark_sra_run_list <- full_df_addit %>%
    filter(
        group != "nuno" & # exclude addit metag provided by collaborators
        group != "bill" & # exclude addit metag provided by collaborators
        group != "metat" & # exclude metat
        group != "sag" & # exclude metat
        depth_group != "epi" | # dark ocean samples from typical regions
        light_avail == "dark" # dark samples from addit seas
    ) %>%
    mutate(dark_sra_run_list = str_c("from_sra/", run_accessions))


write_lines(unique(dark_sra_run_list$dark_sra_run_list), "all_dark_sra_run_v2.txt")

# save final metadata
setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/metadata")
#write_csv(full_df, "metag_metat_sag_v1.csv")

write_csv(full_df_addit, "metag_metat_sag_v2.csv")

