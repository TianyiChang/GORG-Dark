#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

library(tidyverse)
library(RColorBrewer)
library(janitor)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit")

#! noted: baltic, black, ross metag were excluded for frag_recruit analyses due to atypical depth layer
#! noted: NAs in PC values derived from SAGs with few reads mapped to the analyzed metag

#! noted: for a few SAGs (e.g., 'AG-919-M03' and 'AG-908-G06'), they were defined as 'meso' SAGs but embeded in epi sag cloud in the pcoa plot
#! noted: these SAGs have really high relative abundance in some epi samples but absence in others, therefore the difference is not significant.

#======================================#
## Import group info and pcoa outputs ##
#======================================#

#1 grouping info
pref_depth_sag <- 
    read_csv("sag_depth_region_preference/sunlit_dark_major_depth_4_indiv_sags.csv") %>% 
    select(sag, preferred_depth = major_group) %>% 
    mutate(
        preferred_depth = str_replace(preferred_depth, "epi", "EPI"),
        preferred_depth = str_replace(preferred_depth, "meso", "MES"),
        preferred_depth = str_replace(preferred_depth, "bathy", "BAT"),
        preferred_depth = str_replace(preferred_depth, "abysso", "ABY"),
        preferred_depth = str_replace(preferred_depth, "hadal", "HAD"),
        preferred_depth = str_replace_all(preferred_depth, " - ", "-"))

gd_sag_metadata <-
    read_csv("../sag_metadata/v3_SAG_summary_20240320.csv") %>%
    clean_names() %>% 
    mutate(
        depth = ifelse(depth == "omz", 666, depth),
        depth = as.numeric(str_replace(depth, "m|,", "")))

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
relative_eigs <- read_csv("pcoa/sunlit_dark_relative_eigs.csv")
relative_eig1 <- relative_eigs$relative_eig1
relative_eig2 <- relative_eigs$relative_eig2

pc_1_2 <- read_csv("pcoa/sunlit_dark_sag_pc1_2_values.csv") 

#============#
## Join dfs ##
#============#
gd_depth_pc_df <- gd_sag_metadata %>% 
    select(sag) %>% 
    left_join(order_genus_all, by = "sag") %>% 
    left_join(pref_depth_sag, by = "sag") %>% 
    left_join(pc_1_2, by = "sag")

gd_gt_depth_pc_df <- pref_depth_sag %>% 
    left_join(pc_1_2, by = "sag")

# baltic sea, black sea, ross ice shelf gd sags were missing from 'pref_depth_sag',
# add these sags from gd_sag_metadata
gd_gt_depth_pc_df <- gd_depth_pc_df %>% 
    select(sag, preferred_depth, plate, PC_1, PC_2) %>% 
    bind_rows(gd_gt_depth_pc_df) %>% 
    distinct()

# append "depth_niche" column to sag metadata
major_pref_depth <- gd_depth_pc_df %>% 
    count(preferred_depth) %>% 
    mutate(perc = 100 * n /sum(n)) %>% 
    arrange(-perc) %>% 
    filter(n >= 20)

pref_depth_sag_filtered <- pref_depth_sag %>% 
    semi_join(
        major_pref_depth, by = "preferred_depth")

v4_SAG_summary <- read_csv(
    "../sag_metadata/v3_SAG_summary_20240320.csv") %>%
    left_join(pref_depth_sag_filtered, by = c("SAG" = "sag")) %>% 
    mutate(niche_depth = ifelse(
        is.na(preferred_depth), "rare_in_metag", preferred_depth
    )) %>% 
    select(-preferred_depth)

#20241210: update depth_niche after the inclusion of 9 addit hadal metag (n=11)
#20241210: change all Black and Baltic Sea sags depth_niche to NA
#20241218: change depth of AM-164 from 'omz' to 400 as Maria requested
depth_niche_new_df <- v4_SAG_summary %>% 
    select(SAG, niche_depth_new = niche_depth)

v4_SAG_summary_updated <- read_csv("../sag_metadata/v4_SAG_summary_20240625.csv") %>% 
    rename(niche_depth_old = niche_depth) %>% 
    left_join(depth_niche_new_df, by = "SAG") %>% 
    relocate(niche_depth_new, .after = niche_depth_old) %>% 
    mutate(
        niche_depth_new = case_when(
            ocean_province == 'Black Sea' ~ 'NA',
            ocean_province == 'Baltic Sea' ~ 'NA',
            TRUE ~ niche_depth_new
        ),
        depth = ifelse(depth == 'omz', 400, depth),
        depth = as.numeric(str_replace(depth, "m|,", ""))
    )

write_csv(v4_SAG_summary_updated, "../sag_metadata/v5_SAG_summary_20241218.csv")

# append "doubling_hours" column to sag metadata
grodon <- read_tsv("../grodon/gorgd_grodon_combined.tsv") %>% 
    select(SAG = genome, doubling_hours = d)

v4_SAG_summary <- v4_SAG_summary %>% 
    left_join(grodon, by = "SAG")

#! check inconsistent depth layers
v4_SAG_summary <- v4_SAG_summary %>% 
    mutate(
        depth = parse_number(depth),
        sampling_layer = case_when(
            ocean_province == "Black Sea" | ocean_province == "Baltic Sea" | ocean_province == "Ross Ice Shelf" | depth == "omz" | niche_depth == "rare_in_metag" ~ "",
            depth <= 200 ~ "EPI",
            depth > 200 & depth <= 1000 ~ "MES",
            depth > 1000 & depth <= 4000 ~ "BAT",
            depth > 4000 & depth <= 6000 ~ "ABY",
            depth > 6000 ~ "HAD"),
        sampling_layer_niche_depth_mismatch = ifelse(
            str_detect(niche_depth, sampling_layer), "N", "Y"
    ))

write_csv(v4_SAG_summary, "../sag_metadata/v4_SAG_summary_20240625_mismatch.csv")

#! add pa_fl info into the metadata table
pa_fl_df <- read_csv(
    'sag_PA_vs_FL/summary/statistics_w_pa_fl_malaspina_metag_only_1000_4200_sag_depth.csv') %>% 
    select(SAG = sag, pa_fl, pa_fl_pvalue = p.adj)

v4_SAG_summary <- v4_SAG_summary %>% 
    left_join(pa_fl_df, by = "SAG")

write_csv(v4_SAG_summary, "../sag_metadata/v4_SAG_summary_20240823_pafl.csv")

#=====================#
## Generate figures ##
#=====================#

#! noted: check whether the fct_level is correct for 'major_pref_depth'

fct_level = c(
    "EPI", "EPI-MES", "MES", "MES-BAT",
    "BAT", "MES-BAT-ABY", "BAT-ABY",
    "ABY", "ABY-HAD", "HAD")

# raw figures
generate_pcoa_figs <- function(input_df, outfile) {

    # filter out rare combo (e.g., epi-hadal, caused by metag contam)
    major_pref_depth <- {{input_df}} %>% 
        count(preferred_depth) %>% 
        mutate(perc = 100 * n /sum(n)) %>% 
        arrange(-perc) %>% 
        filter(n >= 50)

    filtered_depth_pc_df <- {{input_df}} %>% 
        filter(preferred_depth %in% major_pref_depth$preferred_depth)

    # add factors: region, depth_latitude, preferred depth, open_oceans
    filtered_depth_pc_df <- filtered_depth_pc_df %>%
        mutate(
            preferred_depth = factor(
                preferred_depth,levels = fct_level)
            )

    filtered_depth_pc_df <- mutate(
        filtered_depth_pc_df, preferred_depth = factor(
            preferred_depth, ordered = TRUE)) %>%
        arrange(preferred_depth)

    colour_count <- length(
        unique(filtered_depth_pc_df$preferred_depth))
    #? Spectral (11) or Paired (12) or Bluses (9)
    # colorRampPalette(rev(brewer.pal(11, "Spectral"))) to invert color
    # getPalette = colorRampPalette(brewer.pal(9, "YlGnBu")[3:9])
    getPalette = colorRampPalette(brewer.pal(9, "RdBu"))

    ggplot(filtered_depth_pc_df, aes(x=PC_1, y=PC_2)) +
        geom_point(alpha=0.6, size=1.3,
            aes(colour=factor(preferred_depth))) +
        scale_color_manual(values=
            getPalette(colour_count)) +
        theme_classic() +
        labs(x=str_c("Eig1", " (", relative_eig1, ")"),
            y=str_c("Eig2", " (", relative_eig2, ")")) +
        theme(
            legend.position="bottom",
            legend.text=element_text(size=8, face="bold"),
            legend.title=element_blank()
            ) +
        guides(colour=guide_legend(
            nrow=2,byrow=TRUE,
            override.aes=list(size=4.5)))

    path <- "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/pcoa/raw_figs/"
    ggsave(
        str_c(path, {{outfile}}),
        device="pdf",
        width=6.5, height=6.5)

}

generate_pcoa_figs(
    gd_depth_pc_df,
    "gd_pcoa_preferred_depth_sunlit_dark_metag_rdbu.pdf")

generate_pcoa_figs(
    gd_gt_depth_pc_df,
    "gd_gt_pcoa_preferred_depth_sunlit_dark_metag_rdbu.pdf")

#====================================#
## Upset plot for region preference ##
#====================================#
library(UpSetR)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/sag_depth_region_preference")

region_pref_df <- read_csv("major_region_4_indiv_sags.csv")

express_df <- region_pref_df %>% 
    mutate(
        expression_input = str_replace_all(major_group, " ", "_"),
        expression_input = str_replace_all(expression_input, "_-_", "&")
        ) %>% 
    group_by(expression_input) %>% 
    summarise(count = n()) %>% 
    arrange(-count) %>% 
    filter(!is.na(expression_input))

# express_df_major <- express_df %>% 
#     slice_max(count, n=40)

express_lst <- split(express_df$count, express_df$expression_input)

upset_input <- fromExpression(express_lst)

pdf(
    "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/tmp/sag_region_prefer_upset.pdf",
    width = 9, height = 6)
# the number of bars can be specified with 'nintersects' and 'nsets'
upset(upset_input, 
        nintersects = 50, 
        nsets = 11, 
        order.by = "freq", 
        decreasing = T, 
        mb.ratio = c(0.6, 0.4),
        show.numbers = "no",
        text.scale = 1.1, 
        point.size = 2.8, 
        line.size = 1
      )
dev.off()

