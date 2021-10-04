

import json
with open("gtex_v8.custom.config") as instr:
    gtex_conf = json.load(instr)

configfile: "config.json",

target_tissues = ['as'] #'ebv'] #['as', 'ebv', 'ms', 'wb'] #gtex_conf['tshorts']
target_gxnorm  = ["tpms_qcfilter"]

wildcard_constraints:
    tissue  = "|".join([x for x in target_tissues]),
    gxnorm = "|".join([x for x in target_gxnorm]),
    npeer  = "\d+"

# subworkflow part2:
#     snakefile:
#         "pipe2.snakefile"

tissue_peer_dict = dict()
with open(config['tissue_table']) as ifile:
    for line in ifile:
        if line.startswith("#"): continue
        arr = line.rstrip().split("\t")
        tissue_peer_dict[arr[1]] = arr[3]

target_files = [f"{config['gtex_exprdir']}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.npeer{tissue_peer_dict[tissue]}.txt.gz" for tissue in target_tissues for gxnorm in target_gxnorm]

rule all:
    input:
        target_files

def get_covariate_input(wildcards):
    return f"{config['covariate_dir']}/{gtex_conf[wildcards.tissue]['fullname']}.v8.covariates.txt"

rule prepare_covariates:
    input:
        get_covariate_input
    output:
        outfile=config['covariate_dir']+"/{tissue}.no_peer.txt"
    shell:
        #"grep -iv inferred {input} | awk -f ../scripts/transpose.awk > {output.outfile};"
        "grep -iv inferred {input} > {output.outfile};"

# decompress expr, sort covariates as expr samples, transpose covariate, delete row and col names for covariate
rule prepare_peer:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.txt.gz",
        cov=config['covariate_dir']+"/{tissue}.no_peer.txt"
    output:
        tmpdir   =directory("{workdir}/{gxnorm}/{tissue}_tmp_npeer{npeer}")
    shell:
        "{config[python]} ../scripts/prepare_peer_inputs.py --gx {input.expr} --cov {input.cov} --outdir {output.tmpdir}"

rule calc_peer:
    input:
        tmpdir  = "{workdir}/{gxnorm}/{tissue}_tmp_npeer{npeer}"
    output:
        touch("{workdir}/{gxnorm}/{tissue}_tmp_npeer{npeer}/peer.done")
    shell:
        "peertool -f {input.tmpdir}/gx.tab -n {wildcards.npeer} -c {input.tmpdir}/cov.tab --no_a_out --no_z_out -o {input.tmpdir};"

rule reformat_peer_results:
    input:
        tmpdir  = "{workdir}/{gxnorm}/{tissue}_tmp_npeer{npeer}",
        peerdone= "{workdir}/{gxnorm}/{tissue}_tmp_npeer{npeer}/peer.done",
        origexpr= "{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.txt.gz"
    output:
        outexpr = "{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.npeer{npeer}.txt.gz"
    run:
        import pandas as pd
        orig_df = pd.read_csv(f"{input[2]}", header=0, index_col=0, sep="\t")
        corr_df = pd.read_csv(f"{input[0]}/residuals.csv", header=None, index_col=None, sep=",")
        corr_df.index = list(orig_df.index)
        corr_df.columns = list(orig_df.columns)
        corr_df.to_csv(output[0], header=True, index=True, sep="\t")