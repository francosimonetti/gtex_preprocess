#############################
## Pipeline 2/2 to pre-process GTEx Genotype and RNA-seq expression files
#############################

import json
with open("gtex_v8.custom.config") as instr:
    gtex_conf = json.load(instr)

configfile: "config.json",

target_tissues = gtex_conf['tshorts']
target_gxnorm  = ["tpms_qcfilter", "qn", "tmm", "counts_qcfilter"]

wildcard_constraints:
    tissue  = "|".join([x for x in target_tissues]),
    gxnorm = "|".join([x for x in target_gxnorm])

subworkflow part1:
    snakefile:
        "pipe1.snakefile"


rule all:
    input:
        expand("{mydir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.txt.gz", mydir=config['gtex_exprdir'], tissue=target_tissues, gxnorm=target_gxnorm)

'''
From the raw expression and counts file per tissue, it does some QC filtering. 
By default, it filters genes with <0.1 
'''
rule QC_expression:
    input:
        part1('step1.tpms.done'),
        part1('step1.counts.done'),
    output:
        out="{workdir}/{gxnorm}/{tissue}.{gxnorm}.txt.gz",
    params:
        outdir=config['gtex_exprdir'],
        gexpr_file="{workdir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9.gene_tpm.{tissue}.txt.gz",
        count_file="{workdir}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9.gene_reads.{tissue}.txt.gz"
    shell:
        "{config[python]} ../scripts/gtex_normalization.py --rpkm {params.gexpr_file} --counts {params.count_file} --tissue {wildcards.tissue} --donors {config[donor_file]} --outdir {params.outdir}"

rule filter_genes:
    input:
        infile="{workdir}/{gxnorm}/{tissue}.{gxnorm}.txt.gz"
    output:
        outfile="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.txt.gz"
    params:
        gencode=config['gencode'],
        biotype="protein_coding lncRNA",
    shell:
        "{config[python]} ../scripts/filter_gencode_expr.py --gx {input.infile} --gtf {params.gencode} --biotype {params.biotype} --out {output.outfile}"