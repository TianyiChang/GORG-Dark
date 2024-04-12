#!/home/tchang/miniconda3/envs/shortread/bin/Rscript --vanilla

# this script is used to get normalized abund of SAGs in the collected metag
# sunlit_dark: gorg-dark + gorg-tropics against sunlit and dark metag

library(tidyverse)
library(furrr)
library(Rsamtools)

###################################################
#! obtain number of effective mapped read count ##
###################################################

# noted: plan(multisession) will result in bugs
    
plan(multicore, workers = 12)

get_bam_stats <- function(metag_batch) {

    setwd(str_c(
        "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/",
        metag_batch,
        "/gorg_v3_tropics/aln_filter_bam"
        ))

    get_idxstats <- function(bamfile) {
        idxstatsBam(bamfile) %>%
        mutate(run_accessions = str_replace(
            bamfile, ".bam", ""
        ))
    }

    metag_total_reads <- bind_cols(
        read_csv("../stats/n_mapped_read.csv") %>% 
            select(run_accessions),
        read_csv("../stats/n_total_read.csv") %>% 
            select(n_total_reads)
            ) %>% 
        mutate(run_accessions = str_remove(
            run_accessions, "_wo_Regions.bam$"
        ))

    stats_df <- list.files(pattern = ".bam$", full.names = FALSE) %>%
        future_map_dfr(., get_idxstats) %>%
        mutate(sag = str_replace(seqnames, "_NODE_.*$", "")) %>%
        group_by(run_accessions, sag) %>%
        summarise(mapped_reads = sum(mapped)) %>%
        ungroup() %>% 
        filter(sag != "*") %>%
        left_join(metag_total_reads, by = "run_accessions") # remove results for metag with < 1M reads (the n used for rarefy)

    return(stats_df)

}

combined_stats_df <- bind_rows(
    get_bam_stats("result_4_sra_dark_sunlit"),
    get_bam_stats("result_4_local_dark_sunlit")
    ) %>%
    distinct(sag, run_accessions, .keep_all = TRUE) %>% 
    mutate(sag = str_remove(sag, "_.*"))

##################################
#! obtain taxonomy for each SAG ##
##################################

get_taxonomy <- function(file) {
    setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/gtdbtk")
    read_tsv(file) %>%
    select(user_genome, classification) %>%
    mutate(sag = str_replace(user_genome, "_contigs$", ""))
}

sag_classif <- bind_rows(
    # GORG-Dark
    get_taxonomy("tkv2_gdv3.ar53.summary.tsv"),
    get_taxonomy("tkv2_gdv3.bac120.summary.tsv"),
    # GORG-Tropics
    get_taxonomy("../../simonsproject/gtdbtkv2_gorg_12715/all.ar53.summary.tsv"),
    get_taxonomy("../../simonsproject/gtdbtkv2_gorg_12715/all.bac120.summary.tsv")
    
)

######################################
#! obtain normalized abund per SAG ##
#####################################

# normalizing counts of mapped reads
# the removed sags are not possible surface taxa and lack assembly info

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/sag_metadata")

gd_sag_size <- read_tsv("v3_SAG_summary_20240320.tsv") %>%
    select(sag = SAG, assembly_length)

gt_sag_size <- read_tsv("gorg-tropics_sags_tableS2.tsv") %>% 
    select(sag, assembly_length = assembly_size_bp)

combined_sag_size <- bind_rows(gd_sag_size, gt_sag_size)

stats_df_normalized_count <- combined_sag_size %>% 
    right_join(combined_stats_df, by = "sag") %>%
    mutate(
        mapped_count_perMread_perMbp =
            (1e6 / n_total_reads) *
            (1e6 * mapped_reads / assembly_length)
        ) %>% 
    left_join(sag_classif, by = "sag") %>%
    select(-user_genome) %>%
    # remove the two SAGs with inconsis 16S taxonomy
    filter(
        sag != "AM-259-L06",
        sag != "AM-685-K20"
        )

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/sag_abund_metag")
write_csv(stats_df_normalized_count, "norm_abund_sunlit_dark_sra_collab_gd_gt.csv")
