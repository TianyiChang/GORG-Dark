#!/home/tchang/miniconda3/envs/rgdal/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)
library(rgdal)
library(geosphere)
library(RColorBrewer)

#########################################
## select representative surface metag ##
#########################################

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")
full_df <- read_csv("metadata/metag_metat_sag_v4.csv")

epi_df <- full_df %>%
    filter(depth_group == "epi" &
        (group %in% c("collected", "uncollected")) # include only sra runs
        )

# calculate geographical cluster

xy <- SpatialPointsDataFrame(
    matrix(c(epi_df$longitude, epi_df$latitude), ncol = 2), data.frame(ID = seq_along(1:nrow(epi_df))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

# use the distm function to generate a geodesic distance matrix in meters
mdist <- distm(xy)

# cluster all points using a hierarchical clustering approach
hc <- hclust(as.dist(mdist), method="complete")

# define the distance threshold, d is meter, (277.8 KM = 150 Nautical Miles)
#! 20240401: change to '0', now only derep for longi latitude
d=0

# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
epi_df$clust <- cutree(hc, h=d)

epi_selected <- epi_df %>%
    group_by(clust) %>%
    slice_head(n = 1) %>%
    ungroup() %>% 
    # exclude metag from special region where depth layer cannot be defined
    filter(!ocean_province %in% c("Black Sea", "Baltic Sea", "Ross Ice Shelf"))

# export selected sra_run list
write_lines(epi_selected$run_accessions, "metadata/selected_epi_metag_runs.txt")
write_csv(epi_selected, "metadata/selected_epi_metag.csv")

#! if passing any argument: show the selected on map

if(length(args != 0)) {

    # load shapefile
    setwd("/mnt/scgc/stepanauskas_nfs/simon_lab/tianyi/ocean_map_ggplot")

    wmap <- readOGR(dsn="ne_110m_land", layer="ne_110m_land")
    # convert to dataframe
    wmap_df <- fortify(wmap)

    grat <- readOGR("ne_110m_graticules_all", layer="ne_110m_graticules_15") 
    grat_df <- fortify(grat)

    # create a blank ggplot theme
    theme_opts <- list(theme(panel.grid.minor = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.background = element_blank(),
                            plot.background = element_rect(fill="#e6e8ed"),
                            panel.border = element_blank(),
                            axis.line = element_blank(),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            axis.ticks = element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            plot.title = element_text(size=22)))


    bbox <- readOGR("ne_110m_graticules_all", layer="ne_110m_wgs84_bounding_box") 
    bbox_df<- fortify(bbox)

    # samples<-read.csv("metadata.csv", header = T)
    # samples$Depth<-as.numeric(samples$Depth)
    plot_map <- function(data_frame) {

        ggplot(bbox_df, aes(long,lat, group=group)) + 
        geom_polygon(fill="white") +
        geom_polygon(data=wmap_df, aes(long,lat, group=group, fill=hole)) + 
        #geom_path(data=grat_df, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="lightgrey", size=0.2) +
        geom_point(data=data_frame, aes(longitude, latitude, shape=NULL,
            fill=NULL, group=NULL, size=NULL, color=ocean_province), alpha=I(8/10)) +
        coord_equal() + 
        theme_opts +
        scale_fill_manual(values=c("black", "white"), guide="none")

    }

    setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/fig")

    plot_map(epi_selected)
    ggsave("map_selected_samples_epi.pdf", device="pdf", width=12.5, height=8.25)

} else {
   print("add any args to generate a map")
}