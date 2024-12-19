library(tidyverse)
library(xml2)

args <- commandArgs(trailingOnly = TRUE)

file_info <- file.info(args[1])

# avoid empty file raises error in read_xml
if (file_info$size != 0) {

    xml_data <- read_xml(args[1])

    accession <- xml_data %>% xml_find_all("BioSample") %>% xml_attr("accession")

    metadata_tag <- xml_data %>% xml_find_all("BioSample/Attributes/Attribute") %>% xml_attr('attribute_name')

    metadata_value <- xml_data %>% xml_find_all("BioSample/Attributes/Attribute") %>% xml_text()

    df <- tibble(
            accession = accession,
            metadata_tag = metadata_tag,
            metadata_value = metadata_value) %>% 
        distinct(metadata_tag, .keep_all = TRUE) %>% 
        pivot_wider(
            names_from = metadata_tag,
            values_from = metadata_value)

    write_csv(df, args[2])

}
