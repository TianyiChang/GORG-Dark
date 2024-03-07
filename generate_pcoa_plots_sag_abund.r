#!/home/tchang/miniconda3/envs/shortread/bin/Rscript --vanilla

library(tidyverse)
library(RColorBrewer)
library(janitor)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")

#======================================#
## Import group info and pcoa outputs ##
#======================================#

#1 grouping info
pref_depth_sag <- 
    read_csv("sag_depth_region_preference/major_depth_4_indiv_sags.csv") %>% 
    select(sag, preferred_depth = major_group)

pref_region_sag <- 
    read_csv("sag_depth_region_preference/major_region_4_indiv_sags.csv")

pa_fl_sag <-
    read_csv("sag_PA_vs_FL/summary/statistics_w_pa_fl.csv") %>%
    select(sag, pa_fl) %>%
    distinct()

sag_plate_depth_region <- 
    read.csv("metadata/metag_metat_sag_v3.csv") %>%
    filter(group == "sag") %>%
    as_tibble() %>% 
    select(plate, ocean_province, depth_group)

sag_metadata <-
    read_tsv("../sag_metadata/gdv3_SAG_summary_20231215.tsv") %>%
    clean_names() %>% 
    mutate(
        depth = ifelse(depth == "omz", 666, depth),
        depth = as.numeric(str_replace(depth, "m|,", "")))

depth_latitude <- sag_metadata %>%
    mutate(depth_latitude =
        case_when(
            depth > 1000 & depth <= 4000 ~ "Bathypelagic",
            depth <= 1000 & latitude > 60 ~
                "Mesopelagic & latitude > 60",
            depth <= 1000 & latitude <= 60 ~
                "Mesopelagic & latitude <= 60",
            depth > 4000 & depth <= 6000 ~ "Abyssopelagic",
            depth > 6000 ~ "Hadalpelagic"
            )) %>%
    select(sag, depth_latitude)

get_taxonomy <- function(file) {
    read_tsv(str_c("../gtdbtk/", {{file}})) %>%
    select(user_genome, classification) %>%
    mutate(sag = str_replace(user_genome, "_contigs$", ""))
}

sag_classif <- bind_rows(
    get_taxonomy("tkv2_gdv3.ar53.summary.tsv"),
    get_taxonomy("tkv2_gdv3.bac120.summary.tsv")
    ) %>%
    separate_rows(classification, sep = ";") %>%
    select(-user_genome)

sag_phylum <- sag_classif %>%
    filter(str_detect(classification, "p__"))

sag_orders <- sag_classif %>%
    filter(str_detect(classification, "o__")) %>%
    left_join(sag_phylum, by = "sag") %>%
    mutate(order = str_c(
        classification.x, " (",
        classification.y, ")"
    )) %>%
    select(sag, order)

sag_famlies <- sag_classif %>%
    filter(str_detect(classification, "f__")) %>%
    left_join(sag_orders, by = "sag") %>%
    mutate(
        family = str_c(
            classification, " (", order, ")"),
        family = str_replace(family,
            " \\(p[^\\)]+\\)", "")
    ) %>%
    select(sag, family)

sag_genera <- sag_classif %>%
    filter(str_detect(classification, "g__")) %>%
    left_join(sag_famlies, by = "sag") %>%
    mutate(
        genus = str_c(
            classification, " (", family, ")"),
        genus = str_replace(genus,
            " \\(o[^\\)]+\\)", "")
    ) %>%
    select(sag, genus)

# get major lineage >= 1% and 3%
order_freq <- sag_orders %>%
    count(order) %>%
    arrange(-n) %>%
    mutate(
        prop = 100 * n / sum(n),
        accum_prop = cumsum(prop))
    
order_selected <- order_freq %>%
    filter(prop >= 1) %>%
    left_join(sag_orders, by = "order") %>%
    mutate(prop_cutoff =
        ifelse(prop >= 3, "3per",
        "1per")) %>%
    select(sag, order, prop_cutoff)

family_selected <- sag_famlies %>%
    inner_join(order_selected, by = "sag")

genus_selected <- sag_genera %>%
    inner_join(family_selected, by = "sag")

# taxonomy for all sags
order_genus_all <- sag_genera %>%
    left_join(sag_famlies, by = "sag") %>%
    left_join(sag_orders, by = "sag")

#2 PCoA outputs
relative_eigs <- read_csv("pcoa/relative_eigs.csv")
relative_eig1 <- relative_eigs$relative_eig1
relative_eig2 <- relative_eigs$relative_eig2

pc_1_2 <- read_csv("pcoa/sag_pc1_2_values.csv") %>% 
    select(-plate)

#============#
## Join dfs ##
#============#
combined_df <- sag_metadata %>% 
    select(sag) %>% 
    mutate(plate = str_replace(sag, "(.*)-.*", "\\1")) %>% 
    left_join(sag_plate_depth_region, by = "plate") %>% 
    left_join(depth_latitude, by = "sag") %>% 
    left_join(order_genus_all, by = "sag") %>% 
    left_join(pref_depth_sag, by = "sag") %>% 
    left_join(pc_1_2, by = "sag")

#=====================#
## Generate figures ##
#=====================#

# raw figures
ggplot(combined_df, aes(x=PC_1, y=PC_2)) +
    geom_point(alpha=0.6, size=0.2,
        aes(colour=factor(ocean_province))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=7, face="bold")
        ) +
    guides(colour=guide_legend(
        override.aes=list(size=5)))

ggsave("pcoa/raw_figs/regions.pdf", device="pdf",
    width=6, height=6)

#2
combined_df <- combined_df %>%
    mutate(
        region = 
            case_when(
                ocean_province == "Red Sea" ~ "Red Sea",
                ocean_province == "Mediterranean Sea" ~ "Mediterranean Sea",
                ocean_province == "Arctic Ocean" ~ "Arctic Ocean",
                ocean_province == "Southern Ocean" ~ "Southern Ocean",
                ocean_province == "Southern Ocean" ~ "Southern Ocean",
                ocean_province == "Baltic Sea" ~ "Baltic Sea",
                ocean_province == "Black Sea" ~ "Black Sea",
                TRUE ~ "Open oceans"
            ),
        region = factor(region, levels = c("Open oceans", "Mediterranean Sea", "Red Sea",
            "Southern Ocean", "Arctic Ocean", "Baltic Sea", "Black Sea")),
        depth_latitude = factor(depth_latitude, levels = c(
            "Mesopelagic & latitude <= 60",
            "Mesopelagic & latitude > 60",
            "Bathypelagic", "Abyssopelagic", "Hadalpelagic"
            )),
        preferred_depth = factor(preferred_depth, levels = c(
            "meso", "bathy", "abysso",
            "meso - bathy", "bathy - abysso"))
    )

combined_df_1 <- mutate(combined_df, region = factor(region, ordered = TRUE)) %>%
    arrange(region)

ggplot(combined_df_1, aes(x=PC_1, y=PC_2)) +
    geom_point(alpha=0.7, size=1.2, pch=21,
        aes(colour=region, fill=depth_latitude)) +
    scale_color_manual(values=
        brewer.pal(n = 8, name = "Paired")) +
    scale_fill_manual(values=c("transparent", "gray", "#660000", "black", "red")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6, face="bold"),
        legend.title=element_blank()
        ) +
    guides(
        col=guide_legend(
            nrow=5,
            override.aes=list(size=4.5)),
        fill=guide_legend(
            nrow=5,
            override.aes=list(size=4.5)),
    )

ggsave("pcoa/raw_figs/pcoa_region_depth_latitude.pdf", device="pdf",
    width=6.5, height=6.5)

combined_df_2 <- mutate(combined_df, preferred_depth = factor(preferred_depth, ordered = TRUE)) %>%
    arrange(preferred_depth)

ggplot(combined_df_2, aes(x=PC_1, y=PC_2)) +
    geom_point(alpha=0.7, size=1.2, shape=21,
        aes(fill=preferred_depth, colour=preferred_depth)) +
    scale_fill_manual(values=
        brewer.pal(n = 6, name = "Paired"), na.value = "white") +
    scale_color_manual(values=
        brewer.pal(n = 6, name = "Paired"), na.value = "gray") +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6, face="bold"),
        legend.title=element_blank()
        ) +
    guides(
        fill=guide_legend(
            nrow=5,
            override.aes=list(size=4.5)),
    )

ggsave("pcoa/raw_figs/pcoa_preferred_depth.pdf", device="pdf",
    width=6.5, height=6.5)










#=======================#
## Old version of code ##
#=======================#

possible_surface_sag <-
    read_csv("sag_surface_vs_dark/summary/surface_taxa_statistics.csv") %>%
    select(sag) %>%
    mutate(possible_surface = "Possible surface taxa")

pa_fl_sag <-
    read_csv("sag_PA_vs_FL/summary/statistics_w_pa_fl.csv") %>%
    select(sag, pa_fl) %>%
    distinct()

# multivariate analyses
# points: individ. sags
# df: columns (samples), rows (sags)

normd_abund <-
    read_csv("sag_abund_metag/norm_abund_dark_addit_metag.csv") %>%
    select(sag, run_accessions, mapped_count_perMread_perMbp)

normd_abund_zero <- normd_abund %>%
    group_by(sag) %>%
    summarise(sum = sum(mapped_count_perMread_perMbp)) %>%
    filter(sum == "0")

non_zero_df <- normd_abund %>%
    anti_join(normd_abund_zero, by = "sag") %>%
    pivot_wider(names_from = run_accessions,
        values_from = mapped_count_perMread_perMbp) %>%
    column_to_rownames(var = "sag")

#PCoA from ape
bray_pcoa <- non_zero_df %>%
    vegdist(method = "bray") %>%
    pcoa()

# get eigenvalues for the two major axeses
relative_eig1 <- round(bray_pcoa$values$Relative_eig[1], 2)
relative_eig2 <- round(bray_pcoa$values$Relative_eig[2], 2)

#! get metadata
sag_metadata <-
    read_csv("../GORG-Dark_v1/overview_gorg_dark_v1_20221104.csv") %>%
    mutate(
        depth = ifelse(depth == "omz", 666, depth),
        depth = as.numeric(str_replace(depth, "m", "")))

plate_metadata <-
    read_csv("../sag_metadata/GORG-Dark_SAG_plate_list.csv")

metag_metat_sag <- 
    read.csv("metadata/metag_metat_sag_v2.csv") %>%
    filter(group == "sag") %>% as_tibble()

# suspected grouping    
sag_grouping <- sag_metadata %>%
    mutate(group =
        case_when(
            depth > 2000 ~ "depth >2000m",
            depth <= 2000 & latitude > 60 ~
                "depth <= 2000 & latitude > 60",
            depth <= 2000 & latitude <= 60 ~
                "depth <= 2000 & latitude <= 60"
            )) %>%
    select(sag = sample_id, group)

sag_grouping_2 <- sag_metadata %>%
    mutate(group_2 =
        case_when(
            depth > 1000 & depth <= 4000 ~ "Bathypelagic",
            depth <= 1000 & latitude > 60 ~
                "Mesopelagic & latitude > 60",
            depth <= 1000 & latitude <= 60 ~
                "Mesopelagic & latitude <= 60",
            depth > 4000 & depth <= 6000 ~ "Abyssopelagic",
            depth > 6000 ~ "Hadalpelagic"
            )) %>%
    select(sag = sample_id, group_2)

#! get taxonomy information
get_taxonomy <- function(file) {
    setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/GORG-Dark_v1/gtdb-tk_v2_clean_genomes")
    read_tsv(file) %>%
    select(user_genome, classification) %>%
    mutate(sag = str_replace(user_genome, "_contigs$", ""))
}

sag_classif <- bind_rows(
    get_taxonomy("all.ar53.summary.tsv"),
    get_taxonomy("all.bac120.summary.tsv")
) %>%
separate_rows(classification, sep = ";") %>%
select(-user_genome)

sag_phylum <- sag_classif %>%
    filter(str_detect(classification, "p__"))

sag_orders <- sag_classif %>%
    filter(str_detect(classification, "o__")) %>%
    left_join(sag_phylum, by = "sag") %>%
    mutate(order = str_c(
        classification.x, " (",
        classification.y, ")"
    )) %>%
    select(sag, order)

sag_famlies <- sag_classif %>%
    filter(str_detect(classification, "f__")) %>%
    left_join(sag_orders, by = "sag") %>%
    mutate(
        family = str_c(
            classification, " (", order, ")"),
        family = str_replace(family,
            " \\(p[^\\)]+\\)", "")
    ) %>%
    select(sag, family)

sag_genera <- sag_classif %>%
    filter(str_detect(classification, "g__")) %>%
    left_join(sag_famlies, by = "sag") %>%
    mutate(
        genus = str_c(
            classification, " (", family, ")"),
        genus = str_replace(genus,
            " \\(o[^\\)]+\\)", "")
    ) %>%
    select(sag, genus)

# taxonomy for sags occurrence >= 1%
order_freq <- sag_orders %>%
    count(order) %>%
    arrange(-n) %>%
    mutate(
        prop = 100 * n / sum(n),
        accum_prop = cumsum(prop))
    
order_selected <- order_freq %>%
    filter(prop >= 1) %>%
    left_join(sag_orders, by = "order") %>%
    mutate(prop_cutoff =
        ifelse(prop >= 3, "3per",
        "1per")) %>%
    select(sag, order, prop_cutoff)

family_selected <- sag_famlies %>%
    inner_join(order_selected, by = "sag")

genus_selected <- sag_genera %>%
    inner_join(family_selected, by = "sag")

# taxonomy for all sags
order_genus_all <- sag_genera %>%
    left_join(sag_famlies, by = "sag") %>%
    left_join(sag_orders, by = "sag")

#! get coord positions in the two major axeses
pc_1_2 <- bray_pcoa$vectors[,c(1,2)] %>%
    as_tibble()

df <- normd_abund %>%
    anti_join(normd_abund_zero, by = "sag") %>%
    distinct(sag) %>%
    mutate(
        plate = str_replace(
            sag, "(.*)-.*", "\\1")
            ) %>%
    left_join(possible_surface_sag, by = "sag") %>%
    left_join(pa_fl_sag, by = "sag") %>%
    left_join(metag_metat_sag, by = "plate") %>%
    left_join(sag_grouping, by = "sag") %>%
    left_join(sag_grouping_2, by = "sag") %>%
    left_join(genus_selected, by = "sag")

df_all <- normd_abund %>%
    anti_join(normd_abund_zero, by = "sag") %>%
    distinct(sag) %>%
    left_join(possible_surface_sag, by = "sag") %>%
    left_join(sag_grouping, by = "sag") %>%
    left_join(sag_grouping_2, by = "sag") %>%
    left_join(order_genus_all, by = "sag")

df$possible_surface <-
    replace_na(df$possible_surface, "Possible dark taxa")
df$pa_fl <-    
    replace_na(df$pa_fl, "out of\ndepth range")

df$pa_fl <- factor(df$pa_fl, levels = c("FL",
    "PA", "both", "few_mapped_reads", "out of\ndepth range"))

df_2 <- bind_cols(pc_1_2, df)

df_all$possible_surface <-
    replace_na(df_all$possible_surface, "Possible dark taxa")
df_all_2 <- bind_cols(pc_1_2, df_all)

df_3per <- df_2 %>%
    filter(prop_cutoff == "3per")
df_1per <- df_2 %>%
    filter(prop_cutoff == "1per")

Nitrososphaerales_fam <- df_2 %>%
    filter(str_detect(order, "Nitrososphaerales"))

Nitrosopumilaceae_genu <- df_2 %>%
    filter(str_detect(family, "Nitrosopumilaceae"))

Pelagibacterales_fam <- df_2 %>%
    filter(str_detect(order, "Pelagibacterales"))

Pelagibacteraceae_major_genu_list <- df_2 %>%
    filter(str_detect(family, "Pelagibacteraceae")) %>%
    count(genus) %>%
    mutate(
        count = n,
        prop = 100* count / sum(count)
    ) %>%
    filter(prop >= 0.5)

Pelagibacteraceae_major_genu <- df_2 %>%
    filter(str_detect(family, "Pelagibacteraceae")) %>%
    semi_join(Pelagibacteraceae_major_genu_list, by = "genus")

#! export tables
setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/")

write.csv(non_zero_df, row.names = TRUE,
    "multivariate/normalized_abund_non_zero.csv")
write_csv(df_2, "multivariate/pcoa_coord_variables.csv")

# plot the PCoA result
setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/")

ggplot(df_2, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=0.2,
        aes(colour=factor(possible_surface))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=7, face="bold")
        ) +
    guides(colour=guide_legend(
        override.aes=list(size=5)))

ggsave("multivariate/pcoa_sag_possible_surface.pdf", device="pdf",
    width=6, height=6)

ggplot(df_2, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=0.2,
        aes(colour=factor(group))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=7, face="bold")
        ) +
    guides(colour=guide_legend(nrow=3,byrow=TRUE,
        override.aes=list(size=5)))

ggsave("multivariate/pcoa_sag_depth_latit.pdf", device="pdf",
    width=6, height=6.5)

ggplot(df_2, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=0.2,
        aes(colour=factor(group_2))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=7, face="bold")
        ) +
    guides(colour=guide_legend(nrow=5,byrow=TRUE,
        override.aes=list(size=5)))

ggsave("multivariate/pcoa_sag_depth_latit_2.pdf", device="pdf",
    width=6, height=9)

ggplot(df_1per, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=0.2,
        aes(colour=factor(order))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=5, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=5,byrow=TRUE,
        override.aes=list(size=2.5)))

ggsave("multivariate/pcoa_sag_1per_orders.pdf", device="pdf",
    width=6, height=9)

ggplot(df_3per, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=0.2,
        aes(colour=factor(order))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=7, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=5,byrow=TRUE,
        override.aes=list(size=5)))

ggsave("multivariate/pcoa_sag_3per_orders.pdf", device="pdf",
    width=6, height=9)

ggplot(df_3per, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.7, size=1.2, pch=21,
        aes(colour=order, fill=group)) +
    scale_fill_manual(values=c("white", "red", "black")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    guides(
        col=guide_legend(
            nrow=3,
            override.aes=list(size=4.5)),
        fill=guide_legend(
            nrow=3,
            override.aes=list(size=4.5)),
        ) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6.5, face="bold"),
        legend.title=element_blank()
        )
    
ggsave("multivariate/pcoa_sag_order_depth_latit_2k_depth.pdf", device="pdf",
    width=7.5, height=9)

#! key figure
ggplot(df_3per, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.7, size=1.2, pch=21,
        aes(colour=order, fill=group_2)) +
    scale_fill_manual(values=c("black", "dark red", "white", "dark blue")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    guides(
        col=guide_legend(
            nrow=4,
            override.aes=list(size=4.5)),
        fill=guide_legend(
            nrow=4,
            override.aes=list(size=4.5)),
        ) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6.5, face="bold"),
        legend.title=element_blank()
        )
    
ggsave("multivariate/pcoa_sag_order_depth_latit_1k_depth.pdf", device="pdf",
    width=8, height=9)

#! key figure
ggplot(df_3per, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.7, size=1.2, pch=21,
        aes(colour=order, fill=pa_fl)) +
    scale_fill_manual(values=c("dark red", "black", "dark blue", "grey", "white")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    guides(
        col=guide_legend(
            nrow=3,
            override.aes=list(size=4.5)),
        fill=guide_legend(
            nrow=3,
            override.aes=list(size=4.5)),
        ) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=5.5, face="bold"),
        legend.title=element_blank()
        )

ggsave("multivariate/pcoa_sag_pa_fl.pdf", device="pdf",
    width=7, height=7.5)

# ocean province
colour_count <- length(
    unique(df$ocean_province))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

ggplot(df_2, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=1.3,
        aes(colour=factor(ocean_province))) +
    scale_color_manual(values=
        getPalette(colour_count)) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6.3, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=5,byrow=TRUE,
        override.aes=list(size=4.5)))

ggsave("multivariate/pcoa_sag_ocean_province.pdf", device="pdf",
    width=7, height=7.5)

#! special ocean province
df_niche <- df_2 %>%
    mutate(niche = 
        case_when(
            ocean_province == "Red Sea" ~ "Red Sea",
            ocean_province == "Mediterranean Sea" ~ "Mediterranean Sea",
            ocean_province == "Arctic Ocean" ~ "Arctic Ocean",
            ocean_province == "Southern Ocean" ~ "Southern Ocean",
            ocean_province == "Southern Ocean" ~ "Southern Ocean",
            ocean_province == "Baltic Sea" ~ "Baltic Sea ",
            TRUE ~ "Open oceans"
        ))

ggplot(df_niche, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.3, size=1.3,
        aes(colour=niche)) +
    scale_color_manual(values=
        brewer.pal(n = 7, name = "Set1")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=7, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=2,byrow=TRUE,
        override.aes=list(size=6)))

ggsave("multivariate/pcoa_sag_niche.pdf", device="pdf",
    width=5, height=5.5)

#! key figure, main figure, nich and depth
df_niche <- df_2 %>%
    mutate(niche = 
        case_when(
            ocean_province == "Red Sea" ~ "Red Sea",
            ocean_province == "Mediterranean Sea" ~ "Mediterranean Sea",
            ocean_province == "Arctic Ocean" ~ "Arctic Ocean",
            ocean_province == "Southern Ocean" ~ "Southern Ocean",
            ocean_province == "Southern Ocean" ~ "Southern Ocean",
            ocean_province == "Baltic Sea" ~ "Baltic Sea",
            TRUE ~ "Open oceans"
        ))

df_niche$niche <- factor(
    df_niche$niche, levels = c("Open oceans", "Mediterranean Sea", "Red Sea",
        "Southern Ocean", "Arctic Ocean", "Baltic Sea")
)

df_niche <- mutate(df_niche, niche = factor(niche, ordered = TRUE)) %>%
    arrange(niche)

ggplot(df_niche, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.7, size=1.2, pch=21,
        aes(colour=niche, fill=group_2)) +
    scale_color_manual(values=
        brewer.pal(n = 7, name = "Paired")) +
    scale_fill_manual(values=c("black", "red", "transparent", "grey")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6, face="bold"),
        legend.title=element_blank()
        ) +
    guides(
        col=guide_legend(
            nrow=5,
            override.aes=list(size=4.5)),
        fill=guide_legend(
            nrow=5,
            override.aes=list(size=4.5)),
    )

ggsave("multivariate/pcoa_sag_niche_depth.pdf", device="pdf",
    width=5, height=6.5)


#! 202309: selected
library(patchwork)

plot_niche <- ggplot(df_niche, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.7, size=1, pch=21,
        aes(colour=niche, fill=group_2)) +
    scale_color_manual(values=
        brewer.pal(n = 7, name = "Paired")) +
    scale_fill_manual(values=c("black", "red", "transparent", "grey")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6, face="bold"),
        legend.title=element_blank()
        ) +
    guides(
        col=guide_legend(
            ncol=1,
            override.aes=list(size=4.5)),
        fill=guide_legend(
            ncol=1,
            override.aes=list(size=4.5)),
    )

plot_pa_fl <- ggplot(df_3per, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.7, size=1, pch=21,
        aes(colour=order, fill=pa_fl)) +
    scale_fill_manual(values=c("dark red", "black", "dark blue", "grey", "white")) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    guides(
        col=guide_legend(
            ncol=1,
            override.aes=list(size=4.5)),
        fill=guide_legend(
            ncol=1,
            override.aes=list(size=4.5)),
        ) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6, face="bold"),
        legend.title=element_blank()
        )

plot_niche + plot_pa_fl +
    plot_layout(
        guides = "collect",  # Collect legends
        widths = c(1, 1)    # Adjust relative heights
    ) &
    theme(
        legend.position = "bottom",  # Position legends to the right
        legend.box = "horizontal"  # Arrange legend items horizontally
    )

ggsave("multivariate/pcoa_niche_depth_pa_fl_combined.pdf", device="pdf",
    width=7, height=5)



# major lineages
ggplot(Nitrosopumilaceae_genu, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=1.2,
        aes(colour=factor(genus))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=7, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=5,byrow=TRUE,
        override.aes=list(size=5)))

ggsave("multivariate/pcoa_sag_Nitrosopumilaceae_genu.pdf", device="pdf",
    width=7.5, height=9)

ggplot(Pelagibacteraceae_major_genu, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=1.2,
        aes(colour=factor(genus))) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=7, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=5,byrow=TRUE,
        override.aes=list(size=5)))

ggsave("multivariate/pcoa_sag_Pelagibacteraceae_genu.pdf", device="pdf",
    width=7.5, height=9)

# info about the cluster with the lowest eg2 values
low_eg2_cluster <- df_all_2 %>%
    filter(Axis.1 > -0.08, Axis.1 < -0.015, Axis.2 < -0.43) %>%
    mutate(Plate = str_replace(sag, "-[^-]+$", ""))

low_eg2_cluster_plate <- low_eg2_cluster %>%
    count(Plate) %>%
    left_join(plate_metadata, by = "Plate")

write_csv(low_eg2_cluster, "multivariate/lowest_eg2_cluster_summary.csv")
write_csv(low_eg2_cluster_plate, "multivariate/lowest_eg2_cluster_plate_count.csv")

# other figures requested by RS

low_eg2_cluster_represent <- low_eg2_cluster %>%
    count(genus) %>%
    filter(n >=3, genus != "g__ (f__Pelagibacteraceae)")

df_2_low_eg2_cluster_represent <- df_all_2 %>%
    semi_join(low_eg2_cluster_represent,
        by = "genus")

colour_count <- length(
    unique(df_2_low_eg2_cluster_represent$genus))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

ggplot(df_2_low_eg2_cluster_represent, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=1.3,
        aes(colour=factor(genus))) +
    scale_color_manual(values=
        getPalette(colour_count)) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6.3, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=5,byrow=TRUE,
        override.aes=list(size=4.5)))

ggsave("multivariate/pcoa_sag_major_genu_bottom_cluster.pdf", device="pdf",
    width=7.5, height=9)


# RS: remove cosmopolitan and include cyanob, then re-plot the above
low_eg2_cluster_no_cosmop <- low_eg2_cluster_represent %>%
    select(genus) %>%
    filter(
        genus != "g__ (f__Pelagibacteraceae)",
        genus != "g__Pelagibacter (f__Pelagibacteraceae)",
        genus != "g__Nitromaritima (f__Nitrospinaceae)",
        genus != "g__Nitromaritima (f__Nitrospinaceae)",
        genus != "g__Pseudothioglobus (f__Thioglobaceae)"
        )

cyano_genera <- df_all_2 %>%
    filter(str_detect(order, "Cyano")) %>%
    select(genus)

low_eg2_cluster_wo_cosmop_w_cyano <- bind_rows(
    low_eg2_cluster_no_cosmop,
    cyano_genera
)

df_2_low_eg2_cluster_represent_2 <- df_all_2 %>%
    semi_join(low_eg2_cluster_wo_cosmop_w_cyano,
        by = "genus")

colour_count <- length(
    unique(df_2_low_eg2_cluster_represent_2$genus))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

ggplot(df_2_low_eg2_cluster_represent_2, aes(x=Axis.1, y=Axis.2)) +
    geom_point(alpha=0.6, size=1.3,
        aes(colour=factor(genus))) +
    scale_color_manual(values=
        getPalette(colour_count)) +
    xlim(
        min(df_2_low_eg2_cluster_represent$Axis.1),
        max(df_2_low_eg2_cluster_represent$Axis.1)
        ) +
    ylim(
        min(df_2_low_eg2_cluster_represent$Axis.2),
        max(df_2_low_eg2_cluster_represent$Axis.2)
        ) +
    theme_classic() +
    labs(x=str_c("Eig1", " (", relative_eig1, ")"),
        y=str_c("Eig2", " (", relative_eig2, ")")) +
    theme(
        legend.position="bottom",
        legend.text=element_text(size=6.3, face="bold"),
        legend.title=element_blank()
        ) +
    guides(colour=guide_legend(
        nrow=5,byrow=TRUE,
        override.aes=list(size=4.5)))

ggsave("multivariate/pcoa_sag_major_genu_bottom_cluster_wo_cosmo_w_cyano.pdf", device="pdf",
    width=7.5, height=9)

# # nMDS (too slow for the large dataset)
# set.seed(123) # so the random number generated can be traced

# bray_nmds <- normd_abund_non_zero %>%
#     pivot_wider(names_from = run_accessions,
#         values_from = mapped_count_perMread_perMbp) %>%
#     column_to_rownames(var = "sag") %>%
#     vegdist(method = "bray") %>%
#     metaMDS(k=2, trymax=1000)

# # Adds the variable we later use for coloring to the data frame

# bind_cols(normd_abund_bray$points, normd_abund_df)
