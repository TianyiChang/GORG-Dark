import pandas as pd

run_id_table = pd.read_table('../../sra_run_list.txt')
accession = list(run_id_table.ids)

rule all:
    input:
        expand("../df/{sra}.done", sra=accession)
        
rule retrieve_metadata_xml:
    output:
        "../xml/{sra}.xml"
    conda: "entrez"
    resources:
        mem_mb=20,
        load=50
    threads: 1
    shell:
        '''
        esearch -db sra -query {wildcards.sra}[Accession] \
            | efetch -format xml > {output}
        '''

rule xml2df:
    input:
        xml="../xml/{sra}.xml",
        rscript="../../../../scripts/pcoa_sag_abund/metadata_xml2df.r"
    output:
        touch("../df/{sra}.done")
    params:
        outfile="../df/{sra}.csv"
    conda: "tidyverse"
    resources:
        mem_mb=20,
        load=1
    threads: 1
    shell:
        '''
        Rscript {input.rscript} {input.xml} {params.outfile}
        '''

