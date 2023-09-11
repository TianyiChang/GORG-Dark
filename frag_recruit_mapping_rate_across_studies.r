#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

library(tidyverse)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/")

metadata <- read_csv("metadata/metag_metat_sag_v2.csv") %>%
    mutate(run_accessions = str_replace_all(
        run_accessions, "\\s", ""
    ))

get_n_read <- function(METAG_BATCH) {

    files <- list.files(path = str_c("mapping/", METAG_BATCH),
        pattern = "n_mapped_read.csv|n_total_read.csv",
        full.names = TRUE,
        recursive = TRUE)

    file_list <- list()

    for(i in files) {
        file_list[[i]] <- read_csv(i) %>%
        mutate(
            dataset = str_replace(i,
                "^.*/(.*)/stats.*$", "\\1"),
            metag_batch = METAG_BATCH)
    }

    combined_mapped <- file_list %>%
        bind_rows() %>%
        filter(is.na(n_total_reads)) %>%
        select(-n_total_reads) %>%
        mutate(
            run_accessions = str_replace(run_accessions,
                "_wo_Regions.*", ""),
            run_accessions = str_replace(run_accessions,
                "^(FK21_\\d+)_.*$", "\\1"),
            )

    combined_total <- file_list %>%
        bind_rows() %>%
        filter(!is.na(n_total_reads)) %>%
        select(-n_mapped_reads)

    map_rate_metadata <- combined_mapped %>%
        mutate(mapping_rate = n_mapped_reads / combined_total$n_total_reads) %>%
        filter(dataset != "gorg_v2_concat_old") %>%
        left_join(metadata, by = "run_accessions") %>%
        arrange(desc(group)) %>% select(-group) %>%
        distinct(run_accessions, dataset, .keep_all = TRUE)

    return(map_rate_metadata)

}

combined_map_rate_metadata <- bind_rows(
    get_n_read("dark/w_bwa_aln_filter"),
    get_n_read("addit_collaboraters"),
    get_n_read("addit_seas_particle")
    ) %>%
    filter(
        (is.na(depth_group) | depth_group != "epi") &
        (is.na(light_avail) | light_avail != "light")
        ) %>%
    distinct(run_accessions, dataset, .keep_all = TRUE)

########################
### statistical test ###
########################

# check the class of each variable
lapply(combined_map_rate_metadata, class)

combined_map_rate_metadata$dataset <- factor(combined_map_rate_metadata$dataset)
combined_map_rate_metadata$depth_group <- factor(combined_map_rate_metadata$depth_group)
combined_map_rate_metadata$ocean_province <- factor(combined_map_rate_metadata$ocean_province)

# logit transformation on proportions

# before the transformation: add the 0.01*(the smallest non-zero proportions) to all proportions
col <- combined_map_rate_metadata$mapping_rate

combined_map_rate_metadata$proportion_non_zero <-
    col + (0.01 * min(col[col > 0])) 

# logit transformation: mapping rate
logitTransform <- function(p) { log(p/(1-p)) }

df_logit <- combined_map_rate_metadata %>%
    mutate(logit_rate = logitTransform(proportion_non_zero))

#! part below has been updated
# REML model (REML=TRUE) for better estimating random effects
    # recommended to use on the final selected model 
# ML model for LRT model comparison for estimating fix effects

# LRT for fixed effects
full_model <- lmer(logit_rate ~ dataset +
    depth_group + ocean_province +
    (1 | run_accessions), REML = F, data = df_logit)

model_dataset <- lmer(logit_rate ~ dataset +
    (1 | run_accessions), REML = F, data = df_logit)

model_depth <- lmer(logit_rate ~ depth_group +
    (1 | run_accessions), REML = F, data = df_logit)

model_ocean <- lmer(logit_rate ~ ocean_province +
    (1 | run_accessions), REML = F, data = df_logit)

null_model <- lmer(logit_rate ~ (1 | run_accessions),
    REML = F, data = df_logit)

#! p-value=1.05116e-152
anova(null_model, full_model)
pchisq(758.01, df=14, lower.tail=FALSE)

#! p-value=7.247259e-152
anova(null_model, model_dataset)
pchisq(702.13, df=3, lower.tail=FALSE)

#? p-value=0.0641
anova(null_model, model_depth)
pchisq(7.2586, df=3, lower.tail=FALSE)

#! p-value=0.00283
anova(null_model, model_ocean)
pchisq(23.452, df=8, lower.tail=FALSE)

# emmeans comparison (post-hoc with REML model)
# (https://aosmith.rbind.io/2019/04/15/custom-contrasts-emmeans/)
full_model_reml <- lmer(logit_rate ~ dataset +
    depth_group + ocean_province +
    (1 | run_accessions), REML = T, data = df_logit)

emm_options(opt.digits = FALSE)

emms_dataset <- emmeans(full_model_reml, ~ dataset)
emms_ocean <- emmeans(full_model_reml, ~ ocean_province)

#! between datasets
gorg_v2_concat <- c(1, 0, 0, 0) # assign the position of the vector in 'emms'
outside_acinas_2020 <- c(0, 1, 0, 0)
outside_gorg_tropics <- c(0, 0, 1, 0)
outside_omd_mags <- c(0, 0, 0, 1)

# two-tailed t-test: estimate(-0.3058486); p.value(0.0007419568)
contrast(emms_dataset, method = list(gorg_v2_concat - outside_omd_mags))
2*pt(q=-3.390, df=639, lower.tail=TRUE) #TRUE if t-ratio<0

# two-tailed t-test: estimate(1.83615); p.value(2.262934e-71)
contrast(emms_dataset, method = list(gorg_v2_concat - outside_acinas_2020))
2*pt(q=20.354, df=639, lower.tail=FALSE)

# two-tailed t-test: estimate(2.14195); p.value(7.963649e-90)
contrast(emms_dataset, method = list(gorg_v2_concat - outside_gorg_tropics))
2*pt(q=23.744, df=639, lower.tail=FALSE)

#! between oceans

Arctic <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
Indian <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
Mediterranean <- c(0, 0, 1, 0, 0, 0, 0, 0, 0)
Atlantic_N <- c(0, 0, 0, 1, 0, 0, 0, 0, 0)
Pacific_N <- c(0, 0, 0, 0, 1, 0, 0, 0, 0)
Red_sea <- c(0, 0, 0, 0, 0, 1, 0, 0, 0)
Atlantic_S <- c(0, 0, 0, 0, 0, 0, 1, 0, 0)
Pacific_S <- c(0, 0, 0, 0, 0, 0, 0, 1, 0)
Southern <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)

# two-tailed t-test: estimate(-1.147812); p.value(0.01809293)
contrast(emms_ocean, method = list(Mediterranean - Southern))
2*pt(q=-2.381, df=226.71, lower.tail=TRUE) #TRUE if t-ratio<0

#################
### box plots ###
#################

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/")

ggplot(data=df_logit,
    aes(x=dataset, y=logit_rate, fill=dataset)) +
    geom_boxplot() +
    geom_jitter(width=0.25, alpha=1, size=0.3) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("between_datasets/v2/box_transf_dataset.pdf",
    device="pdf", height=5, width=5)

ggplot(data=df_logit,
    aes(x=ocean_province, y=logit_rate, fill=dataset)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("between_datasets/v2/box_transf_ocean_by_dataset.pdf",
    device="pdf", height=6, width=8)

ggplot(data=df_logit,
    aes(x=depth_group, y=logit_rate, fill=dataset)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("between_datasets/v2/box_transf_depth_by_dataset.pdf",
    device="pdf", height=6, width=8)

# original mapping rate
ggplot(data=df_logit,
    aes(x=dataset, y=mapping_rate, fill=dataset)) +
    geom_boxplot() +
    geom_jitter(width=0.25, alpha=1, size=0.3) +
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("between_datasets/v2/box_orig_dataset.pdf",
    device="pdf", height=5, width=5)

#! update 202309
setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/")
df <- read_csv("tables/mapping_rate_all_metag_across_studies.csv")

df_modified <- df %>%
    mutate(
        region_group = ifelse(
            (ocean_province == "Baltic Sea" |
            ocean_province == "Black Sea" |
            ocean_province == "Arctic Ocean" |
            ocean_province == "Southern Ocean" |
            ocean_province == "Red Sea" |
            ocean_province == "Mediterranean Sea" |
            ocean_province == "Ross Ice Shelf"),
            ocean_province, "Open Ocean"
        )
    ) %>%
    filter(dataset != "GORG-Tropics")

df_modified$region_group <- factor(df_modified$region_group,
    levels = c("Black Sea", "Baltic Sea", "Red Sea",
    "Mediterranean Sea", "Ross Ice Shelf", "Southern Ocean",
    "Arctic Ocean", "Open Ocean")
)

df_modified$dataset <- factor(df_modified$dataset,
    levels = c("GORG-Dark", "OMD-MAGs", "MDeep-MAGs"))

ggplot(data=df_modified,
    aes(x=region_group, y=mapping_rate, fill=dataset)) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values=c("#FF5733", "#FFC300", "#1ABC9C")) +
    theme_classic() +
    labs(y = "Fragment recruitment rate", x= "") +
    theme(
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 8),
        legend.position = c(0.5, 0.95),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1)
        ) +
    guides(fill = guide_legend(nrow = 1))

ggsave("fig/box_mapping_rate_dataset.pdf",
    device="pdf", height=5, width=5)

