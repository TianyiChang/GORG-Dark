#!/home/tchang/miniconda3/envs/shortread/bin/Rscript --vanilla

# this script is used to get normalized abund of SAGs in the collected metag

library(tidyverse)
library(furrr)
library(Rsamtools)

###################################################
#! obtain number of effective mapped read count ##
###################################################

# noted: plan(multisession) will result in bugs
    
plan(multicore, workers = 15)

get_bam_stats <- function(metag_batch) {

    setwd(str_c(
        "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/",
        metag_batch,
        "/gorg_v3_concat/aln_filter_bam"
        ))

    get_idxstats <- function(bamfile) {
        idxstatsBam(bamfile) %>%
        mutate(run_accessions = str_replace(
            bamfile, ".bam", ""
        ))
    }

    filter_list <- read_csv("../post_qc_summary/metag_post_qc_read.csv") %>%
        filter(survived_reads < 1000000) #! changes needed, currently throw away hadal metag

    stats_df <- list.files(pattern = ".bam$", full.names = FALSE) %>%
        future_map_dfr(., get_idxstats) %>%
        mutate(sag = str_replace(seqnames, "_NODE_.*$", "")) %>%
        group_by(run_accessions, sag) %>%
        summarise(mapped_reads = sum(mapped)) %>%
        filter(sag != "*") %>%
        anti_join(filter_list, by = "run_accessions") # remove results for metag with < 1M reads (the n used for rarefy)

    return(stats_df)

}

combined_stats_df <- bind_rows(
    get_bam_stats("result_4_sra_metag"),
    get_bam_stats("result_4_local_metag")
    ) %>%
    distinct(sag, run_accessions, .keep_all = TRUE)

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
    get_taxonomy("tkv2_gdv3.ar53.summary.tsv"),
    get_taxonomy("tkv2_gdv3.bac120.summary.tsv")
)

######################################
#! obtain normalized abund per SAG ##
#####################################

# normalizing counts of mapped reads
# the removed sags are not possible surface taxa and lack assembly info

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/sag_metadata")

stats_df_normalized_count <- read_tsv("gdv3_SAG_summary_20231215.tsv") %>%
    select(sag = SAG, assembly_length) %>%
    right_join(combined_stats_df, by = "sag") %>%
    mutate(mapped_count_perMread_perMbp = 
        1e6 * mapped_reads / assembly_length) %>%
    left_join(sag_classif, by = "sag") %>%
    select(-user_genome) %>%
    filter(sag != "AH-888-F19",
        sag != "AH-988-D20",
        sag != "AM-262-O02")

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/sag_abund_metag")
write_csv(stats_df_normalized_count, "norm_abund_dark_collab_metag_gdv3.csv")