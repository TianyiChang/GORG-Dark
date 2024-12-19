#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

# this script generate a world map showing read mapping rates
# for metag libraries as indiv pie charts

library(tidyverse)
library(rsvg)
library(ggimage)
library(maps)
library(scatterpie)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")


# Get clust info for analyzed metag

metadata <- read_csv("metadata/metag_metat_sag_v4.csv") %>%
    mutate(run_accessions = str_replace_all(run_accessions,
        "\\s", ""))

metadata_dark_metag <- metadata %>%
    # filter for metag metad
    filter(group != "sag" & group != "metat") %>% 
    # filter for dark ocean metag: conditions for specialized envs (e.g., Black Sea) and others
    filter(light_avail == "dark" | depth_group != "epi") %>%
    arrange(
        -metat_same_sample, -metat_same_clust,
        -surf_same_clust, group, .by_group = TRUE) %>%
    distinct(run_accessions, .keep_all = TRUE)

sag_metadata <- metadata %>%
    filter(group == "sag") %>%
    arrange(
        -metat_same_sample, -metat_same_clust,
        -surf_same_clust, group, .by_group = TRUE)


get_mean_mapping_rate_by_clust <- function(PATH, depth_group) {
    
    n_mapped_reads <- read_csv(
            str_c(PATH, "gorg_v4_concat/stats/n_mapped_read.csv")) %>%
        mutate(
            run_accessions = str_replace(run_accessions,
            "_wo_Regions.*", ""),
            run_accessions = str_replace(run_accessions,
            "^(FK21_\\d+)_.*$", "\\1"))

    n_total_reads <- read_csv(
        str_c(PATH, "gorg_v4_concat/stats/n_total_read.csv")) %>%
        mutate(run_accessions = n_mapped_reads$run_accessions)
        
    map_rate <- n_total_reads %>%
        left_join(n_mapped_reads, by = "run_accessions") %>%
        mutate(
            mapping_rate = n_mapped_reads / n_total_reads,
            n_unmapped_reads = n_total_reads - n_mapped_reads) %>%
        left_join(metadata_dark_metag, by = "run_accessions") %>%
        filter(
            !is.na(latitude) &
            (is.na(light_avail) | light_avail == "dark"))

    #! note: Assign all Black Sea, Baltic Sea, and Ross metag to 'meso'
    if({{depth_group}} == "meso") {
        mean_mapping_rate_by_clust <- map_rate %>%
            filter(depth_group == "meso" | (light_avail == "dark" & is.na(depth_group))) %>% # calculate aver mapping rate by location ('clust') and depth group
            group_by(clust) %>%
            mutate(mean_m_rate = mean(mapping_rate)) %>%
            ungroup() %>%
            distinct(clust, .keep_all = TRUE) %>%
            mutate(
                n_mapped_reads = mapping_rate,
                n_unmapped_reads = 1 - mapping_rate)
    } else {
        mean_mapping_rate_by_clust <- map_rate %>%
            filter(depth_group == {{depth_group}}) %>% # calculate aver mapping rate by location ('clust') and depth group
            group_by(clust) %>%
            mutate(mean_m_rate = mean(mapping_rate)) %>%
            ungroup() %>%
            distinct(clust, .keep_all = TRUE) %>%
            mutate(
                n_mapped_reads = mapping_rate,
                n_unmapped_reads = 1 - mapping_rate)
    }

    
    return(mean_mapping_rate_by_clust)

}


sra_dark_meso <- get_mean_mapping_rate_by_clust("result_4_sra_dark/", "meso")
sra_dark_bathy <- get_mean_mapping_rate_by_clust("result_4_sra_dark/", "bathy")
sra_dark_abysso <- get_mean_mapping_rate_by_clust("result_4_sra_dark/", "abysso")
sra_dark_hadal <- get_mean_mapping_rate_by_clust("result_4_sra_dark/", "hadal")

local_dark_meso <- get_mean_mapping_rate_by_clust("result_4_local_dark/", "meso")
local_dark_bathy <- get_mean_mapping_rate_by_clust("result_4_local_dark/", "bathy")
local_dark_abysso <- get_mean_mapping_rate_by_clust("result_4_local_dark/", "abysso")
local_dark_hadal <- get_mean_mapping_rate_by_clust("result_4_local_dark/", "hadal")

sag_metadata_by_location <- sag_metadata %>%
    distinct(longitude, latitude, .keep_all =  TRUE)

# plot world ocean map
# read shapefile
create_figure_color_point <- function(df1, df2, outfile) {

    combined_meso <- bind_rows(df1, df2)

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
        p + geom_point(
                aes(x=longitude, y=latitude, color = mapping_rate),
                size=3, data=combined_meso) +
                scale_color_continuous(trans="reverse")

    ggsave({{outfile}}, device = "pdf", width = 12.5, height = 8.25)

}


create_figure_color_point(sra_dark_meso, local_dark_meso, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_meso_color_pts.pdf")
create_figure_color_point(sra_dark_bathy, local_dark_bathy, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_bathy_color_pts.pdf")
create_figure_color_point(sra_dark_abysso, local_dark_abysso, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_abysso_color_pts.pdf")
create_figure_color_point(sra_dark_hadal, local_dark_hadal, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_hadal_color_pts.pdf")


# create_pie_figure <- function(df1, df2, outfile) {

#     combined_meso <- bind_rows(df1, df2)

#         world_map <- map_data('world') %>%
#             rename(Longitude = long, Latitude = lat) %>%
#             ggplot(aes(Longitude, Latitude)) +
#             geom_map(map=map_data('world'), aes(map_id=region), fill="light grey", color="light grey") +
#             coord_quickmap() +
#             scale_x_continuous(breaks=seq(-180, 180, 60)) +
#             scale_y_continuous(breaks=seq(-90, 90, 30)) +
#             theme_classic() +
#             theme(
#                 axis.title=element_text(size=18,color='black',face='bold'),
#                 axis.text=element_text(size=12,color='black',face='bold')
#             )

#         # add icons representing sampling location
#         sag_metadata_by_location$icon <- "boat"

#         being_transparent <- function(img) {
#             magick::image_fx(img, expression = "0.5*a", channel = "alpha")
#         }

#         p <- world_map + geom_icon(data=sag_metadata_by_location,
#             aes(x=longitude, y=latitude, image=icon),
#                 size=0.02, asp=2.5, image_fun=being_transparent)

#         # add pies
#         p + geom_scatterpie(
#                 aes(x=longitude, y=latitude, r=3),
#                 data=combined_meso,
#                 cols=c("n_mapped_reads", "n_unmapped_reads"),
#                 color="black", alpha=0.75
#                 ) +
#                 scale_fill_manual(values=c("red", "white"), guide="none")

#     ggsave({{outfile}}, device = "pdf", width = 12.5, height = 8.25)

# }

# create_pie_figure(sra_dark_meso, local_dark_meso, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_meso.pdf")
# create_pie_figure(sra_dark_bathy, local_dark_bathy, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_bathy.pdf")
# create_pie_figure(sra_dark_abysso, local_dark_abysso, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_abysso.pdf")
# create_pie_figure(sra_dark_hadal, local_dark_hadal, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig/glob_mean_mapping_rate_hadal.pdf")


