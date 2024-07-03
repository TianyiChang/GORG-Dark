#!/home/tchang/miniconda3/envs/tidyverse/bin/Rscript --vanilla

library(tidyverse)
library(rstatix)

# library(lme4)
# library(lmerTest)
# library(pbkrtest)
# library(emmeans)

setwd("/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/")

metadata <- read_csv("metadata/metag_metat_sag_v4.csv") %>%
    mutate(run_accessions = str_replace_all(
        run_accessions, "\\s", ""
    ))

get_n_read <- function(path) {

    files <- list.files(path = {{path}},
        pattern = "n_mapped_read.csv|n_total_read.csv",
        full.names = TRUE,
        recursive = TRUE)

    file_list <- list()

    for(i in files) {
        file_list[[i]] <- read_csv(i) %>%
        mutate(
            dataset = str_replace(i,
                "^.*/(.*)/stats.*$", "\\1"))
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
        mutate(
            mapping_rate = n_mapped_reads / combined_total$n_total_reads) %>%
        left_join(metadata, by = "run_accessions") %>%
        arrange(desc(group)) %>%
        select(-group) %>%
        distinct(
            run_accessions, dataset, .keep_all = TRUE) %>% 
    
    return(map_rate_metadata)

}

#todo: un-comment the filter for result_4_local_dark and result_4_sra_dark
combined_map_rate_metadata <- bind_rows(
    get_n_read("result_4_local_dark"),
    get_n_read("result_4_sra_dark")
    ) %>%
    filter(
        (is.na(depth_group) | depth_group != "epi") &
        (is.na(light_avail) | light_avail != "light")
        ) %>%
    distinct(run_accessions, dataset, .keep_all = TRUE) %>% 
    mutate(
        depth_group = str_replace(depth_group, "epi", "EPI"),
        depth_group = str_replace(depth_group, "meso", "MES"),
        depth_group = str_replace(depth_group, "bathy", "BAT"),
        depth_group = str_replace(depth_group, "abysso", "ABY"),
        depth_group = str_replace(depth_group, "hadal", "HAD"),
        depth_group = ifelse(is.na(depth_group), "Special", depth_group),
        depth_group = factor(depth_group, levels = c("EPI", "MES", "BAT", "ABY", "HAD", "Special")))
    

########################
### statistical test ###
########################

# # check the class of each variable
# lapply(combined_map_rate_metadata, class)

# combined_map_rate_metadata$dataset <- factor(combined_map_rate_metadata$dataset)
# combined_map_rate_metadata$depth_group <- factor(combined_map_rate_metadata$depth_group)
# combined_map_rate_metadata$ocean_province <- factor(combined_map_rate_metadata$ocean_province)

# # logit transformation on proportions

# # before the transformation: add the 0.01*(the smallest non-zero proportions) to all proportions
# col <- combined_map_rate_metadata$mapping_rate

# combined_map_rate_metadata$proportion_non_zero <-
#     col + (0.01 * min(col[col > 0])) 

# # logit transformation: mapping rate
# logitTransform <- function(p) { log(p/(1-p)) }

# df_logit <- combined_map_rate_metadata %>%
#     mutate(logit_rate = logitTransform(proportion_non_zero))

# #! part below has been updated
# # REML model (REML=TRUE) for better estimating random effects
#     # recommended to use on the final selected model 
# # ML model for LRT model comparison for estimating fix effects

# # LRT for fixed effects
# full_model <- lmer(logit_rate ~ dataset +
#     depth_group + ocean_province +
#     (1 | run_accessions), REML = F, data = df_logit)

# model_dataset <- lmer(logit_rate ~ dataset +
#     (1 | run_accessions), REML = F, data = df_logit)

# model_depth <- lmer(logit_rate ~ depth_group +
#     (1 | run_accessions), REML = F, data = df_logit)

# model_ocean <- lmer(logit_rate ~ ocean_province +
#     (1 | run_accessions), REML = F, data = df_logit)

# null_model <- lmer(logit_rate ~ (1 | run_accessions),
#     REML = F, data = df_logit)

# #! p-value=1.05116e-152
# anova(null_model, full_model)
# pchisq(758.01, df=14, lower.tail=FALSE)

# #! p-value=7.247259e-152
# anova(null_model, model_dataset)
# pchisq(702.13, df=3, lower.tail=FALSE)

# #? p-value=0.0641
# anova(null_model, model_depth)
# pchisq(7.2586, df=3, lower.tail=FALSE)

# #! p-value=0.00283
# anova(null_model, model_ocean)
# pchisq(23.452, df=8, lower.tail=FALSE)

# # emmeans comparison (post-hoc with REML model)
# # (https://aosmith.rbind.io/2019/04/15/custom-contrasts-emmeans/)
# full_model_reml <- lmer(logit_rate ~ dataset +
#     depth_group + ocean_province +
#     (1 | run_accessions), REML = T, data = df_logit)

# emm_options(opt.digits = FALSE)

# emms_dataset <- emmeans(full_model_reml, ~ dataset)
# emms_ocean <- emmeans(full_model_reml, ~ ocean_province)

# #! between datasets
# gorg_v2_concat <- c(1, 0, 0, 0) # assign the position of the vector in 'emms'
# outside_acinas_2020 <- c(0, 1, 0, 0)
# outside_gorg_tropics <- c(0, 0, 1, 0)
# outside_omd_mags <- c(0, 0, 0, 1)

# # two-tailed t-test: estimate(-0.3058486); p.value(0.0007419568)
# contrast(emms_dataset, method = list(gorg_v2_concat - outside_omd_mags))
# 2*pt(q=-3.390, df=639, lower.tail=TRUE) #TRUE if t-ratio<0

# # two-tailed t-test: estimate(1.83615); p.value(2.262934e-71)
# contrast(emms_dataset, method = list(gorg_v2_concat - outside_acinas_2020))
# 2*pt(q=20.354, df=639, lower.tail=FALSE)

# # two-tailed t-test: estimate(2.14195); p.value(7.963649e-90)
# contrast(emms_dataset, method = list(gorg_v2_concat - outside_gorg_tropics))
# 2*pt(q=23.744, df=639, lower.tail=FALSE)

# #! between oceans

# Arctic <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
# Indian <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
# Mediterranean <- c(0, 0, 1, 0, 0, 0, 0, 0, 0)
# Atlantic_N <- c(0, 0, 0, 1, 0, 0, 0, 0, 0)
# Pacific_N <- c(0, 0, 0, 0, 1, 0, 0, 0, 0)
# Red_sea <- c(0, 0, 0, 0, 0, 1, 0, 0, 0)
# Atlantic_S <- c(0, 0, 0, 0, 0, 0, 1, 0, 0)
# Pacific_S <- c(0, 0, 0, 0, 0, 0, 0, 1, 0)
# Southern <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)

# # two-tailed t-test: estimate(-1.147812); p.value(0.01809293)
# contrast(emms_ocean, method = list(Mediterranean - Southern))
# 2*pt(q=-2.381, df=226.71, lower.tail=TRUE) #TRUE if t-ratio<0

# logit transformation: mapping rate
combined_map_rate_metadata$dataset <- factor(combined_map_rate_metadata$dataset)
combined_map_rate_metadata$depth_group <- factor(combined_map_rate_metadata$depth_group)
combined_map_rate_metadata$ocean_province <- factor(combined_map_rate_metadata$ocean_province)

# logit transformation on proportions

# before the transformation: add the 0.01*(the smallest non-zero proportions) to all proportions
col <- combined_map_rate_metadata$mapping_rate

combined_map_rate_metadata$proportion_non_zero <-
    col + (0.01 * min(col[col > 0])) 

logitTransform <- function(p) { log(p/(1-p)) }

df_logit <- combined_map_rate_metadata %>%
    mutate(logit_rate = logitTransform(proportion_non_zero))

# filter out groups without enough duplications
df_logit_no_duplic_depth_group <- df_logit %>% 
    group_by(depth_group) %>% 
    count(dataset) %>% 
    filter(n < 5)

df_logit_no_duplic_ocean_province <- df_logit %>% 
    group_by(ocean_province) %>% 
    count(dataset) %>% 
    filter(n < 5)

# games-howell test
games_gt_gd_sunlit_dark_frag_recru_depth <- df_logit %>%
    anti_join(df_logit_no_duplic_depth_group, by = "depth_group") %>% 
    filter(!is.na(depth_group)) %>% # exclude ross, baltic, and black
    group_by(depth_group) %>% 
    games_howell_test(
        logit_rate ~ dataset, conf.level = 0.95, detailed = TRUE)

games_gt_gd_sunlit_dark_frag_recru_region <- df_logit %>%
    anti_join(df_logit_no_duplic_ocean_province, by = "ocean_province") %>% 
    filter(!is.na(ocean_province)) %>% # exclude ross, baltic, and black
    group_by(ocean_province) %>% 
    games_howell_test(
        logit_rate ~ dataset, conf.level = 0.95, detailed = TRUE)

write.csv(
    games_gt_gd_sunlit_dark_frag_recru_depth,
    "between_datasets/v4/games_gt_gd_sunlit_dark_frag_recru_depth.csv")

write.csv(
    games_gt_gd_sunlit_dark_frag_recru_region,
    "between_datasets/v4/games_gt_gd_sunlit_dark_frag_recru_region.csv")

#? Tukey test
tukey_gt_gd_sunlit_dark_frag_recru_depth <- df_logit %>%
    anti_join(df_logit_no_duplic_depth_group, by = "depth_group") %>% 
    filter(!is.na(depth_group)) %>% # exclude ross, baltic, and black
    group_by(depth_group) %>% 
    tukey_hsd(
        logit_rate ~ dataset, conf.level = 0.95, detailed = TRUE)

tukey_gt_gd_sunlit_dark_frag_recru_region <- df_logit %>%
    anti_join(df_logit_no_duplic_ocean_province, by = "ocean_province") %>% 
    filter(!is.na(ocean_province)) %>% # exclude ross, baltic, and black
    group_by(ocean_province) %>% 
    tukey_hsd(
        logit_rate ~ dataset, conf.level = 0.95, detailed = TRUE)

write.csv(
    tukey_gt_gd_sunlit_dark_frag_recru_depth,
    "between_datasets/v4/tukey_gt_gd_sunlit_dark_frag_recru_depth.csv")

write.csv(
    tukey_gt_gd_sunlit_dark_frag_recru_region,
    "between_datasets/v4/tukey_gt_gd_sunlit_dark_frag_recru_region.csv")

#################
### box plots ###
#################

# ggplot(data=df_logit,
#     aes(x=dataset, y=logit_rate, fill=dataset)) +
#     geom_boxplot() +
#     geom_jitter(width=0.25, alpha=1, size=0.3) +
#     theme_classic() +
#     theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

# ggsave("between_datasets/v2/box_transf_dataset.pdf",
#     device="pdf", height=5, width=5)

# ggplot(data=df_logit,
#     aes(x=ocean_province, y=logit_rate, fill=dataset)) +
#     geom_boxplot() +
#     theme_classic() +
#     theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

# ggsave("between_datasets/v2/box_transf_ocean_by_dataset.pdf",
#     device="pdf", height=6, width=8)

# ggplot(data=df_logit,
#     aes(x=depth_group, y=logit_rate, fill=dataset)) +
#     geom_boxplot() +
#     theme_classic() +
#     theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

# ggsave("between_datasets/v2/box_transf_depth_by_dataset.pdf",
#     device="pdf", height=6, width=8)


# original mapping rate
combined_map_rate_metadata$dataset <- factor(
    combined_map_rate_metadata$dataset, levels = c(
        "gorg_v4_concat", "outside_omd_mags", "gorg_v4_omd",
        "outside_acinas_2020", "outside_gorg_tropics"))


combined_map_rate_metadata %>% 
    ggplot(aes(x=dataset, y=mapping_rate, fill=dataset)) +
        geom_boxplot(outlier.shape=NA) +
        xlab("") +
        theme_classic() +
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("between_datasets/v4/box_orig_dataset.pdf",
    device="pdf", height=5, width=5)

combined_map_rate_metadata %>% 
    anti_join(df_logit_no_duplic_ocean_province, by = "ocean_province") %>% 
    ggplot(aes(x=ocean_province, y=mapping_rate, fill=dataset)) +
        geom_boxplot(outlier.shape=NA) +
        xlab("") +
        theme_classic() +
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("between_datasets/v4/box_orig_ocean_province.pdf",
    device="pdf", height=5, width=8)

combined_map_rate_metadata %>% 
    anti_join(df_logit_no_duplic_depth_group, by = "depth_group") %>% 
    ggplot(aes(x=depth_group, y=mapping_rate, fill=dataset)) +
        geom_boxplot(outlier.shape=NA) +
        xlab("") +
        theme_classic() +
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

ggsave("between_datasets/v4/box_orig_depth_group.pdf",
    device="pdf", height=5, width=8)
