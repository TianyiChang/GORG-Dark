#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

# this script generate a world map showing read mapping rates
# for metag libraries as indiv pie charts

library(tidyverse)
library(rsvg)
library(ggimage)
library(maps)
library(scatterpie)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")

metadata <- read_csv("metadata/metag_metat_sag_v2.csv") %>%
    mutate(run_accessions = str_replace_all(run_accessions,
        "\\s", ""))

metadata_2 <- metadata %>%
    filter(group != "sag" & group != "metat" &
        (is.na(depth_group) | depth_group != "epi")) %>%
    arrange(-metat_same_sample, -metat_same_clust,
        -surf_same_clust, group, .by_group = TRUE) %>%
    distinct(run_accessions, .keep_all = TRUE)

sag_metadata <- metadata %>%
    filter(group == "sag") %>%
    arrange(-metat_same_sample, -metat_same_clust,
        -surf_same_clust, group, .by_group = TRUE)

get_mean_mapping_rate_by_clust <- function(PATH) {
    
    n_mapped_reads <-
        read_csv(str_c
            ("mapping/", PATH, "/stats/n_mapped_read.csv")) %>%
        mutate(
            run_accessions = str_replace(run_accessions,
            "_wo_Regions.*", ""),
            run_accessions = str_replace(run_accessions,
            "^(FK21_\\d+)_.*$", "\\1"),
            )

    n_total_reads <- read_csv(str_c
        ("mapping/", PATH, "/stats/n_total_read.csv")) %>%
        mutate(run_accessions = n_mapped_reads$run_accessions)
        
    map_rate <- n_total_reads %>%
        left_join(n_mapped_reads, by = "run_accessions") %>%
        mutate(
            mapping_rate = n_mapped_reads / n_total_reads,
            n_unmapped_reads = n_total_reads - n_mapped_reads
            ) %>%
        left_join(metadata_2, by = "run_accessions") %>%
        filter(!is.na(latitude) &
            (is.na(light_avail) | light_avail == "dark"))

    mean_mapping_rate_by_clust <- map_rate %>%
        group_by(clust) %>%
        mutate(mean_m_rate = mean(mapping_rate)) %>%
        ungroup() %>%
        distinct(clust, .keep_all = TRUE) %>%
        mutate(
            n_mapped_reads = mapping_rate,
            n_unmapped_reads = 1 - mapping_rate
        )
    
    return(mean_mapping_rate_by_clust)

}

dark <- get_mean_mapping_rate_by_clust("dark/w_bwa_aln_filter/gorg_v2_concat")

# remove acinas particle metag that already included in 'dark'
addit_seas <- get_mean_mapping_rate_by_clust("addit_seas_particle/gorg_v2_concat") %>%
    filter(group != "collected", group != "uncollected")

addit_collab <- get_mean_mapping_rate_by_clust("addit_collaboraters/gorg_v2_concat")

sag_metadata_by_location <- sag_metadata %>%
    distinct(longitude, latitude, .keep_all =  TRUE)

combined_df <- bind_rows(
    dark, addit_collab, addit_seas
)

# plot world ocean map
# read shapefile

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig")

pdf("mean_mapping_rate_dark_addit_metag.pdf", width=12.5, height=8.25)

    world_map <- map_data('world') %>%
        rename(Longitude = long, Latitude = lat) %>%
        ggplot(aes(Longitude, Latitude)) +
        geom_map(map=map_data('world'), aes(map_id=region), fill="light grey", color="light grey") +
        coord_quickmap() +
        scale_x_continuous(breaks=seq(-180, 180, 60)) +
        scale_y_continuous(breaks=seq(-90, 90, 30)) +
        theme_classic() +
        theme(
            axis.title=element_text(size=18,color='black',face='bold'),
            axis.text=element_text(size=12,color='black',face='bold')
        )

    # add icons representing sampling location
    sag_metadata_by_location$icon <- "boat"

    being_transparent <- function(img) {
        magick::image_fx(img, expression = "0.5*a", channel = "alpha")
    }

    p <- world_map + geom_icon(data=sag_metadata_by_location,
        aes(x=longitude, y=latitude, image=icon),
            size=0.02, asp=2.5, image_fun=being_transparent)

    # add pies
    p + geom_scatterpie(
            aes(x=longitude, y=latitude, r=3),
            data=combined_df,
            cols=c("n_mapped_reads", "n_unmapped_reads"),
            color="black", alpha=0.75
            ) +
            scale_fill_manual(values=c("red", "white"), guide="none")

dev.off()