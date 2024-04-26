library(tidyverse)
library(data.table)
library(janitor)

get_col_by_keys <- function(path2file) {

    fread({{path2file}}) %>% 
        clean_names() %>% 
        select(
            contains("accession") |
            contains("oxygen") | contains("o2") |
            contains("nitrogen") | contains("nitrate") | contains("nitrite") | contains("ammonium") | 
            contains("sulfur") | contains("sulfate") | contains("sulfite") |
            contains("temp") | contains("salinity") |
            contains("depth") | contains("region") |
            contains("size fraction") | contains("size_fraction") |
            contains("longitude") | contains("latitude") |
            contains("lat_lon") | contains("lon_lat")
            ) %>% 
        mutate_all(as.character)

}


combined_df <- list.files(
    path = "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/metadata/metag_metadata_240420/df",
    pattern = ".csv",
    full.names = TRUE
) %>% 
map_dfr(get_col_by_keys)

# further filter the df to simplify the process
filtered_df <- combined_df %>% 
    select(
        accession, diss_oxygen, oxygen_sensor,
        temp, temperature, salinity_sensor,
        nitrate_sensor, ammonium, nitrate, nitrite
    )

write_csv(filtered_df, "/mnt/scgc/stepanauskas_nfs/projects/gorg-dark/frag_recruit/metadata/metag_metadata_240420/key_env_metad_temp.csv")

nrow(filtered_df)

filtered_df %>% 
    filter(
        !is.na(diss_oxygen) |
        !is.na(oxygen_sensor) 
    )