#!/home/tchang/miniconda3/envs/shortread/bin/Rscript --vanilla

library(tidyverse)
library(furrr)
library(Rsamtools)
library(rstatix)
library(scales)

###################################################
#! obtain number of effective mapped read count ##
###################################################

# noted: plan(multisession) will result in bugs
#! noted: particle-associated metag have already been included in 'mapping/dark'
#! therefore, no need to be re-imported from 'mapping/addit_seas_particle'
    
plan(multicore)

setwd(str_c(
    "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/mapping/",
    "dark/w_bwa_aln_filter/gorg_v2_concat/aln_filter_bam"
    ))

get_idxstats <- function(bamfile) {
    idxstatsBam(bamfile) %>%
    mutate(run_accessions = str_replace(
        bamfile, ".bam", ""
    ))
}

filter_list <- read_csv("../post_qc_summary/metag_post_qc_read.csv") %>%
    filter(survived_reads < 1000000)

stats_df <- list.files(pattern = ".bam$", full.names = FALSE) %>%
    future_map_dfr(., get_idxstats) %>%
    mutate(sag = str_replace(seqnames, "_NODE_.*$", "")) %>%
    group_by(run_accessions, sag) %>%
    summarise(mapped_sag = sum(mapped)) %>%
    filter(sag != "*") %>%
    anti_join(filter_list, by = "run_accessions") # remove results for metag with < 1M reads (the n used for rarefy)

get_taxonomy <- function(file) {
    setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/GORG-Dark_v1/gtdb-tk_v2_clean_genomes")
    read_tsv(file) %>%
    select(user_genome, classification) %>%
    mutate(sag = str_replace(user_genome, "_contigs$", ""))
}

sag_classif <- bind_rows(
    get_taxonomy("all.ar53.summary.tsv"),
    get_taxonomy("all.bac120.summary.tsv")
)

######################################################
#! statistics for possible particle-associated taxa ##
######################################################

# normalizing counts of mapped reads
# the removed sags are not possible surface taxa and lack assembly info


setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/GORG-Dark_v1")

stats_df_normalized_count <- read_csv("overview_gorg_dark_v1_20221104.csv") %>%
    select(sag = sample_id, assembly_length) %>%
    right_join(stats_df, by = "sag") %>%
    mutate(mapped_count_RPKM = 
        1e3 * mapped_sag / assembly_length) %>%
    left_join(sag_classif, by = "sag") %>%
    select(-user_genome) %>%
    filter(sag != "AH-888-F19",
        sag != "AH-988-D20",
        sag != "AM-262-O02")

# only keep the est. abund for metag in depth range bathy-abysso,
# same as the range for PA metag
#! MAYBE SHOULD ONLY USE metag FROM 'The Malaspina Expedition 2010', since the size-fraction is clear

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/metadata")

metag_metadata <- read_csv("metag_metat_sag_v1.csv")

malaspina_norm_abund <- read_csv("metag_to_collect/acinas_malaspina_particle.csv") %>%
    mutate(
        pa_fl = ifelse(`filter size (μm)` == 0.8,
            "PA", "FL")
    ) %>%
    select(
        run_accessions = ENA_Run_Accession_ID,
        pa_fl
        ) %>%
    select(run_accessions, pa_fl) %>%
    left_join(stats_df_normalized_count, by = "run_accessions") %>%
    filter(!is.na(sag))

malaspina_norm_abund_per_sag <- malaspina_norm_abund %>%
    group_by(sag) %>%
    summarise(median_mapped_count_RPKM =
        median(mapped_count_RPKM))

sag_metadata <-
    read_csv("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark//GORG-Dark_v1/overview_gorg_dark_v1_20221104.csv") %>%
    mutate(
        depth = ifelse(depth == "omz", 666, depth),
        depth = as.numeric(str_replace(depth, "m", ""))) %>%
    select(sag = sample_id, depth)


# statistical test

no_read_mapped_bathy_abysso <- malaspina_norm_abund %>%
    group_by(sag) %>%
    summarise(n=sum(mapped_count_RPKM)) %>%
    ungroup() %>%
    filter(n == "0")

malaspina_norm_abund_no_zeros <-
    malaspina_norm_abund %>%
    anti_join(no_read_mapped_bathy_abysso, by = "sag") %>%
    mutate(
        sag = as.factor(sag),
        pa_fl = as.factor(pa_fl),
        mapped_n_dummy =
            mapped_count_RPKM +
            rnorm(length(mapped_count_RPKM), 0.05, 0.0001)
        )

games_howell_summary <- malaspina_norm_abund_no_zeros %>%
    group_by(sag) %>%
    games_howell_test(mapped_count_RPKM ~ pa_fl, conf.level = 0.95, detailed = TRUE) %>%
    left_join(sag_classif, by = "sag") %>%
    left_join(malaspina_norm_abund_per_sag, by = "sag") %>%
    left_join(sag_metadata, by = "sag")

#! noted: estimate = group2-group1

median_abund_all_sags <-
    median(games_howell_summary$median_mapped_count_RPKM)

games_howell_flag <- games_howell_summary %>%
    arrange(-estimate, p.adj) %>%
    mutate(
        order = str_replace(classification,
            "^(.*o__.+);f__.*$", "\\1"),
        pa_fl = case_when(
            p.adj < 0.05 & estimate > 0 ~ "PA",
            p.adj < 0.05 & estimate < 0 ~ "FL",
            p.adj >= 0.05 &
                median_mapped_count_RPKM >= median_abund_all_sags
                    ~ "PA+FL",
            TRUE ~ "Few mapped reads"
        )
    ) %>%
    select(sag, classification,
        median_mapped_count_RPKM,
        est.difference = estimate, p.adj, p.adj.signif,
        sample_size_FL = n1, sample_size_PA = n2,
        order, pa_fl, depth)

# select only SAGs recovered from the same depth range as PA metag if the SAGs were defined as 'free-living'
games_howell_flag_same_depth_range <- games_howell_flag %>%
    filter(pa_fl == "PA" | (depth > 1000 & depth < 4200)) 

prevalent_pa_orders <- games_howell_flag %>%
    filter(pa_fl == "PA") %>%
    count(order) %>%
    arrange(-n)

highly_pa_enriched_orders <- games_howell_flag %>%
    filter(pa_fl == "PA") %>%
    group_by(order) %>%
    summarise(
        est_diff_mean = mean(est.difference)
    ) %>%
    arrange(-est_diff_mean)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/sag_PA_vs_FL/summary")

write_csv(games_howell_flag, "statistics_w_pa_fl_malaspina_metag_only_all_sag_depth_range.csv")
write_csv(games_howell_flag_same_depth_range, "statistics_w_pa_fl_malaspina_metag_only_1000_4200_sag_depth.csv")
write_csv(games_howell_summary, "all_sags_statistics_malaspina_metag_only.csv")
write_csv(prevalent_pa_orders, "prevalent_pa_orders_malaspina_metag_only.csv")
write_csv(highly_pa_enriched_orders, "highly_pa_enriched_orders_malaspina_metag_only.csv")
write_csv(malaspina_norm_abund_no_zeros, "bathy_abysso_pa_fl_metag_norm_abund_malaspina_metag_only.csv")


#################################################
#! boxplots for possible particle-attaced taxa ##
#################################################

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/sag_PA_vs_FL")

malaspina_norm_abund_no_zeros <- read_csv("summary/bathy_abysso_pa_fl_metag_norm_abund_malaspina_metag_only.csv")
games_howell_flag <- read_csv("summary/statistics_w_pa_fl_malaspina_metag_only.csv")

pa_sar11_norm_abund <- games_howell_flag %>%
    filter(
        pa_fl == "PA",
        str_detect(classification, "o__Pelagibacterales")
    ) %>%
    select(sag, pa_fl_sag = pa_fl) %>%
    left_join(malaspina_norm_abund_no_zeros,
        by = "sag") %>%
    select(sag, metag = run_accessions, pa_fl_metag = pa_fl,
        mapped_count_RPKM, classification) %>%
    mutate(
        genus = str_replace(classification,
            ".*;g__(.*);s__.*", "\\1")
    )

setwd("./figure")

p <- ggplot(pa_sar11_norm_abund, aes(x = genus, y = mapped_count_RPKM)) +
        geom_boxplot(
            aes(fill = pa_fl_metag),
                position = position_dodge(0.9),
                outlier.shape = NA 
                ) +
    ylim(c(0, 0.35)) +
    xlab(NULL) + ylab("Reads Per Kilobase Million") +
    scale_fill_manual(values = c("#999999", "#E69F00"),
        name ="", labels =
            c("0.22-0.8 µm", "0.8-20 µm")) +
    theme_classic()

ggsave(p, file="sar11_pa_genera_box_malaspina_metag_only.pdf",
    device="pdf", width = 6, height = 6)


subc2_highly_pa_enriched <- games_howell_flag %>%
    filter(str_detect(
        classification, "AG-414-E02")) %>%
    arrange(-est.difference, p.adj) %>%
    slice_head(n=20) %>%
    mutate(fct_order = c(1:20)) %>%
    select(sag, fct_order) %>%
    left_join(pa_sar11_norm_abund, by = "sag")

subc2_highly_pa_enriched$sag <- fct_reorder(subc2_highly_pa_enriched$sag,
        subc2_highly_pa_enriched$fct_order, min)

subc2_sign_pa_enriched <- games_howell_flag %>%
    filter(str_detect(
        classification, "AG-414-E02")) %>%
    arrange(p.adj) %>%
    slice_head(n=20) %>%
    mutate(fct_order = c(1:20)) %>%
    select(sag, fct_order) %>%
    left_join(pa_sar11_norm_abund, by = "sag")

subc2_sign_pa_enriched$sag <- fct_reorder(subc2_sign_pa_enriched$sag,
        subc2_sign_pa_enriched$fct_order, min)

p2 <- ggplot(subc2_highly_pa_enriched, aes(x = sag, y = mapped_count_RPKM)) +
        geom_boxplot(
            aes(fill = pa_fl_metag),
                position = position_dodge(0.9),
                outlier.size = 0.5
                ) +
    xlab(NULL) + ylab("Reads Per Kilobase Million") +
    scale_fill_manual(values = c("#999999", "#E69F00"),
        name ="", labels =
            c("0.22-0.8 µm", "0.8-20 µm")) +
    theme_classic() +
    theme(
        legend.position = c(0.85, 0.9),
        axis.text.x = element_text(face="bold", angle = 45,
            color="dark blue",size=6, vjust=1, hjust = 1)
            ) 

ggsave(p2, file="sar11_subc2_high_pa_enriched_box_malaspina_metag_only.pdf",
    device="pdf", width = 6.5, height = 4)

p3 <- ggplot(subc2_sign_pa_enriched, aes(x = sag, y = mapped_count_RPKM)) +
        geom_boxplot(
            aes(fill = pa_fl_metag),
                position = position_dodge(0.9),
                outlier.size = 0.5
                ) +
    xlab(NULL) + ylab("Reads Per Kilobase Million") +
    scale_fill_manual(values = c("#999999", "#E69F00"),
        name ="", labels =
            c("0.22-0.8 µm", "0.8-20 µm")) +
    theme_classic() +
    theme(
        legend.position = c(0.85, 0.9),
        axis.text.x = element_text(face="bold", angle = 45,
            color="dark blue",size=6, vjust=1, hjust = 1)
            ) 

ggsave(p3, file="sar11_subc2_sign_pa_enriched_box_malaspina_metag_only.pdf",
    device="pdf", width = 6.5, height = 4)


# surface_taxa_df <- idxstats_df %>%
#     semi_join(surface_taxa_statistics, by = "sag") %>%
#     left_join(sag_classif, by = "sag")

# surface_sag <- unique(surface_taxa_df$sag)
# surface_sag_plots <- list()

# for(sag_ in surface_sag) {

#     data <- filter(surface_taxa_df, sag == sag_)
#     classification <- unique(data$classification)
#     classification <- str_replace_all(classification, ";", "; ")

#     surface_sag_plots[[sag_]] <-
#         ggplot(data, aes(x=mapped_sag, fill=layer, color=layer)) +
#             geom_histogram(aes(y=2*after_stat(density)/sum(after_stat(density))),
#                 position="identity", alpha=0.3) +
#             scale_color_manual(values = c("#999999", "#E69F00")) +
#             scale_fill_manual(values = c("#999999", "#E69F00")) +
#             scale_y_continuous(labels=percent_format()) +
#             xlab("Number of mapped reads") +
#             ylab("Percentage of libraries") +
#             theme_classic() +
#             labs(title = sag_,
#                 subtitle = str_wrap(classification, width = 60)) +
#             theme(plot.subtitle=element_text(size=8, hjust=0.5,
#                 face="bold", colour="maroon"))

#     ggsave(surface_sag_plots[[sag_]], file=str_c("hist_", sag_, ".pdf"),
#         device="pdf", width = 5, height = 5)

# }
