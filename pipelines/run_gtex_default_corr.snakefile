
import json
with open("gtex_v8.custom.config") as instr:
    gtex_conf = json.load(instr)

configfile: "config.json",
target_tissues = ['as'] #['as', 'ebv', 'ms', 'wb'] #gtex_conf['tshorts']

wildcard_constraints:
    tissue  = "|".join([x for x in target_tissues])

rule all:
    input:
        # expand(config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.reformat.txt", tissue=target_tissues)
        expand(config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.covariate_corrected.txt", tissue=target_tissues)

def get_covariate_input(wildcards):
    return f"{config['covariate_dir']}/{gtex_conf[wildcards.tissue]['fullname']}.v8.covariates.txt"

def get_expr_input(wildcards):
    return f"{config['gtexportal_dir']}/{gtex_conf[wildcards.tissue]['fullname']}.v8.normalized_expression.bed.gz"


rule reformat_expression:
    input:
        get_expr_input
        #config['gtexportal_dir']+f"/{tissue_name}.v8.normalized_expression.bed.gz"
        #lambda wc: get_expression_name(wc, append=".v8.normalized_expression.bed.gz")
    output:
        config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.reformat.txt"
    shell:
        "zcat {input} | cut -f 4- > {output}"

rule cclm_correction:
    input:
        expr=config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.reformat.txt",
        cov=get_covariate_input
    output:
        config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.covariate_corrected.txt"
    shell:
        "Rscript ../scripts/cclm.R {input.expr} {input.cov} {output}"

rule gzip_files:
    input:
        expr_reformat=config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.reformat.txt",
        expr_corr=config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.covariate_corrected.txt"
    output:
        expr_reformatgz=config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.reformat.txt.gz",
        expr_corrgz=config['gtexportal_dir']+"/{tissue}.v8.normalized_expression.covariate_corrected.txt.gz"
    shell:
        "gzip {input.expr_reformat}; gzip {input.expr_corr}"