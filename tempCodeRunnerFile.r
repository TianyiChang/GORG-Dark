filtered_df %>% 
    filter(all(is.na(oxygen_sensor), is.na(diss_oxygen)))