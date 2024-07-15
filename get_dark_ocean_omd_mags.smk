import pandas as pd
import re
import glob

files = glob.glob('../../omd_mags/*.fa')
biosamples = [re.sub(r'^../../omd_mags/[^_]+_([^_]+)_.*$', '\\1', biosample) for biosample in files]

rule all:
    input:
        expand("../df/{biosample}.done", biosample=biosamples)
        
rule retrieve_metadata_xml:
    output:
        "../xml/{biosample}.xml"
    conda: "entrez"
    resources:
        mem_mb=20,
        load=50
    threads: 1
    shell:
        '''
        esearch -db BioSample -query {wildcards.biosample}[Accession] \
            | efetch -format xml > {output}
        '''

rule xml2df:
    input:
        xml="../xml/{biosample}.xml",
        rscript="../../../../../scripts/pcoa_sag_abund/metadata_xml2df.r"
    output:
        touch("../df/{biosample}.done")
    params:
        outfile="../df/{biosample}.csv"
    conda: "tidyverse"
    resources:
        mem_mb=20,
        load=1
    threads: 1
    shell:
        '''
        Rscript {input.rscript} {input.xml} {params.outfile}
        '''

