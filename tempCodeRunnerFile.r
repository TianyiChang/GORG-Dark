test <- sag_norm_abund_ori %>% 
    mutate(run_accessions = str_remove(
        run_accessions, "_S\\d{2}_HB27filtered")
        )
        
filter(test, is.na(mapped_count_perMread_perMbp))
