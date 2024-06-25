#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

# this script generate taxa composition barcharts based on gtdb-tk and depth-niche specified using statistical test

library(tidyverse)
library(RColorBrewer)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark")

#===============#
# Import tables #
#===============#

gd_metad <- read_csv("sag_metadata/v5_SAG_summary_20240523.csv")
gt_gd_depth_niche <- read_csv("frag_recruit/sag_depth_region_preference/sunlit_dark_major_depth_4_indiv_sags.csv")

gt_gd_depth_niche <- gt_gd_depth_niche %>% 
    mutate(
        major_group = str_replace(major_group, "epi", "EPI"),
        major_group = str_replace(major_group, "meso", "MES"),
        major_group = str_replace(major_group, "bathy", "BAT"),
        major_group = str_replace(major_group, "abysso", "ABY"),
        major_group = str_replace(major_group, "hadal", "HAD"),
        major_group = str_replace_all(major_group, " - ", "-"))

get_taxonomy <- function(path2file) {
    read_tsv({{path2file}}) %>%
    select(user_genome, classification) %>%
    mutate(sag = str_replace(user_genome, "_contigs$", ""))
}

sag_classif <- bind_rows(
    get_taxonomy("gtdbtk/tkv2_gdv3.ar53.summary.tsv"),
    get_taxonomy("gtdbtk/tkv2_gdv3.bac120.summary.tsv"),
    get_taxonomy("../simonsproject/gtdbtkv2_gorg_12715/all.bac120.summary.tsv"),
    get_taxonomy("../simonsproject/gtdbtkv2_gorg_12715/all.ar53.summary.tsv")
    ) %>%
    separate_rows(classification, sep = ";") %>%
    select(-user_genome)

sag_phylum <- sag_classif %>%
    filter(str_detect(classification, "p__")) %>% 
    rename(phylum = classification) %>% 
    mutate(phylum = str_remove(phylum, "p__"))

sag_order <- sag_classif %>%
    filter(str_detect(classification, "o__")) %>% 
    rename(order = classification) %>% 
    mutate(order = str_remove(order, "o__"))

special_region <- gd_metad %>% 
    filter(
        ocean_province == "Black Sea" |
        ocean_province == "Baltic Sea" |
        ocean_province == "Ross Ice Shelf") %>% 
    select(
        major_group = ocean_province,
        classification = gtdb_classification) %>% 
    separate_rows(classification, sep = ";")

special_region_phylum <- special_region %>% 
    filter(str_detect(classification, "p__")) %>% 
    rename(phylum = classification) %>% 
    mutate(phylum = str_remove(phylum, "p__"))
    
special_region_order <- special_region %>% 
    filter(str_detect(classification, "o__")) %>% 
    rename(order = classification) %>% 
    mutate(order = str_remove(order, "o__"))
    
#==================#
# Generate figures #
#==================#

generate_plot <- function(df1, df2, type1, type2, Type, cutoff) {

    taxa_niche <- gt_gd_depth_niche %>% 
        left_join({df1}, by = "sag") %>% 
        filter(!is.na({{type1}})) %>% #! SAGs wo gtdktk results are excluded
        group_by(major_group) %>% 
        mutate(niche_count = n()) %>% 
        ungroup() %>% 
        # filter(major_group == "epi" | major_group == "meso" |major_group == "bathy" |major_group == "hadal" |major_group == "abysso" |major_group == "hadal") %>% #! RS suggested to only keep single layers
        filter(niche_count > 50) %>% #! SAGs assigned to rare niche type are excluded
        bind_rows({df2})

    relat_abund <- taxa_niche %>% 
        count(major_group, {{type1}}) %>% 
        group_by(major_group) %>% 
        mutate(relat_abund = n / sum(n))  %>% 
        ungroup()

    major_taxa <- relat_abund %>% 
        filter(relat_abund >= {{cutoff}}) %>% 
        pull({{type1}}) %>% 
        unique()

    fct_level = c(
        "EPI", "EPI-MES", "MES", "MES-BAT",
        "BAT", "MES-BAT-ABY", "BAT-ABY",
        "ABY", "ABY-HAD", "HAD",
         "Ross Ice Shelf", "Black Sea", "Baltic Sea")

    # fct_level = c(
    #     "EPI",  "MES", "BAT", "ABY", "HAD",
    #      "Ross Ice Shelf", "Black Sea", "Baltic Sea") #! RS suggested to only keep single layers

    df <- relat_abund %>% 
        mutate(taxa = ifelse(
            {{type1}} %in% major_taxa, {{type1}}, "Others")) %>% 
        group_by(major_group, taxa) %>% 
        summarise(relat_abund = sum(relat_abund)) %>% 
        arrange(-relat_abund) %>% 
        ungroup() %>% 
        mutate(
        taxa = fct_inorder(taxa),
        taxa = fct_relevel(taxa, "Others", after = Inf),
        major_group = factor(major_group, levels = fct_level))

    colourCount <- length(unique(df$taxa))

    ggplot(df, aes(fill=taxa, y=relat_abund, x=major_group)) + 
        geom_bar(position="fill", stat="identity", colour='grey40', linewidth = 0.25) +
        scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(colourCount)) +
        labs(x = "", y = "Relative abudance") +
        labs(fill = str_c({{Type}}, " (GTDB)")) +
        theme_classic() + 
        theme(
            legend.key.size = unit(0.12, 'inch'), legend.text=element_text(size=5), legend.title=element_text(size=8,face="bold"),
            axis.text.y = element_text(size=6), axis.title=element_text(size=8,face="bold"),
            axis.line = element_line(colour = "black"), axis.text.x = element_text(face="bold", angle = 20, color="dark blue",size=6, vjust=1, hjust = 1))

    ggsave(str_c(
        "frag_recruit/taxa_comp/taxa_comp_depth_niche_", {{type2}}, ".pdf"),
        device = "pdf", height = 3.5, width = 4)

}

generate_plot(sag_phylum, special_region_phylum, phylum, "phylum", "Phylum", 0.03)
generate_plot(sag_order, special_region_order, order, "order", "Order", 0.04)
