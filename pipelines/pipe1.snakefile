#############################
## Pipeline 1/2 to pre-process GTEx Genotype and RNA-seq expression files
#############################

import json
with open("gtex_v8.custom.config") as instr:
    gtex_conf = json.load(instr)

configfile: "config.json",

rule all:
    input:
        'step1.tpms.done',
        'step1.counts.done',

'''
Splits the GTEx Portal raw expression files for every tissue listed in the pheno_file
It requires a table with all these tissues and acronyms, one per row, in tissue_table
'''
rule split_expr_tpms:
    input:
        config['expression_gct']
    output:
        touch('step1.tpms.done'),
    params:
        outdir=directory(config["gtex_exprdir"]),
        tissuefile = config['tissue_table'],
        phenofile  = config['pheno_file'],
        donorfile  = config['donor_file']
    shell:
        "{config[python]} ../scripts/split_expr_matrix_by_tissue.py --input {input} --outdir {params.outdir} --tissuefile {params.tissuefile} --pheno {params.phenofile} --donors {params.donorfile}"

rule split_expr_reads:
    input:
        config['counts_gct']
    output:
        touch('step1.counts.done'),
    params:
        outdir=directory(config["gtex_exprdir"]),
        tissuefile = config['tissue_table'],
        phenofile  = config['pheno_file'],
        donorfile  = config['donor_file']
    shell:
        "{config[python]} ../scripts/split_expr_matrix_by_tissue.py --input {input} --outdir {params.outdir} --tissuefile {params.tissuefile} --pheno {params.phenofile} --donors {params.donorfile}"

