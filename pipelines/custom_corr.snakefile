

import json
with open("gtex_v8.custom.config") as instr:
    gtex_conf = json.load(instr)

configfile: "config.json",

target_tissues = ['ebv'] #['as', 'ebv', 'ms', 'wb'] #gtex_conf['tshorts']
target_gxnorm  = ["tpms_qcfilter"]

wildcard_constraints:
    tissue  = "|".join([x for x in target_tissues]),
    gxnorm = "|".join([x for x in target_gxnorm]),
    npeer  = "\d+",
    pca    = "\d+"

# subworkflow part2:
#     snakefile:
#         "pipe2.snakefile"

# PCA   = [3, 3]
# NPEER = [3, 0]
# KNN   = [0, 30]

PCA   = [0,  3,  3]
NPEER = [0,  0,  10]
KNN   = [30, 30, 30]

# targets = expand([config['gtex_exprdir']+"/{{gxnorm}}/{{tissue}}.{{gxnorm}}.pc_lncrna.pca{pca}_npeer{npeer}_knn{k}.txt.gz"], zip, pca=PCA, npeer=NPEER, k=KNN)
# print(targets)


##### INCOMPLETE
### Chequear:
#   - Correr peer acá y peer en mi compu de casa con el script PEER.R . Tendría que dar muy similar. No va a dar igual por el proceso iterativo de PEER


rule all:
    input:
        #expand(targets, gxnorm=target_gxnorm, tissue=target_tissues )
        #expand(config['covariate_dir']+"/{tissue}.pca0.txt", tissue=target_tissues)
        expand(config['gtex_exprdir']+"/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}.txt.gz", pca =PCA, npeer=NPEER, gxnorm=target_gxnorm, tissue=target_tissues),
        #config['gtex_exprdir']+"/{{gxnorm}}/{{tissue}}.{{gxnorm}}.pc_lncrna.pca{pca}_npeer{npeer}_knn{k}.txt.gz"

def get_covariate_input(wildcards):
    return f"{config['covariate_dir']}/{gtex_conf[wildcards.tissue]['fullname']}.v8.covariates.txt"

rule prepare_pca_covariates:
    input:
        get_covariate_input
    output:
        outfile=config['covariate_dir']+"/{tissue}.pca{pca}.txt"
    shell:
        "grep -iv inferred {input} > {output.outfile}.tmp ; head -n $(( {wildcards.pca} + 1 )) {output.outfile}.tmp > {output.outfile}; rm {output.outfile}.tmp"
        #"grep -iv inferred {input} | awk -f ../scripts/transpose.awk > {output.outfile};"

rule cclm_correction:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.txt.gz",
        cov=config['covariate_dir']+"/{tissue}.pca{pca}.txt"
    output:
        "{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}.txt.gz"
    shell:
        "Rscript ../scripts/cclm.R {input.expr} {input.cov} {output}"

# decompress expr, sort covariates as expr samples, transpose covariate, delete row and col names for covariate
rule prepare_peer:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}.txt.gz",
        cov=config['covariate_dir']+"/{tissue}.pca{pca}.txt"
    output:
        tmpdir=directory("{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}")
    shell:
        "{config[python]} ../scripts/prepare_peer_inputs.py --gx {input.expr} --cov {input.cov} --outdir {output.tmpdir}"

rule calc_peer:
    input:
        tmpdir = "{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}"
    output:
        touch("{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}/peer.done")
    shell:
        "if [ {wildcards.npeer} -ne 0 ] ; then peertool -f {input.tmpdir}/gx.tab -n {wildcards.npeer} -c {input.tmpdir}/cov.tab --no_a_out --no_z_out --e_pb 10 --e_pa 0.1 --a_pb 0.01 --a_pa 0.001 -o {input.tmpdir};"
        "else sed 's/\t/,/g' {input.tmpdir}/gx.tab > {input.tmpdir}/residuals_t.csv; awk -f ../scripts/transpose_csv.awk {input.tmpdir}/residuals_t.csv >  {input.tmpdir}/residuals.csv; fi;"

# rule calc_peer_nocov:
#     input:
#         tmpdir = "{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}"
#     output:
#         touch("{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}/peer.done")
#     shell:
#         ### PEER with input covariates or with previous cclm of those covariates is the same! yay!
#         "if [ {wildcards.npeer} -ne 0 ] ; then peertool -f {input.tmpdir}/gx.tab -n {wildcards.npeer} --no_a_out --no_z_out --e_pb 10 --e_pa 0.1 --a_pb 0.01 --a_pa 0.001 -o {input.tmpdir};"
#         "else sed 's/\t/,/g' {input.tmpdir}/gx.tab > {input.tmpdir}/residuals_t.csv; awk -f ../scripts/transpose_csv.awk {input.tmpdir}/residuals_t.csv >  {input.tmpdir}/residuals.csv; fi;"

## peertool -f gx.tab -n 3 --no_a_out --no_z_out -c cov.tab -o test --e_pb 10 --e_pa 0.1 --a_pb 0.01 --a_pa 0.001

rule reformat_peer_results:
    input:
        tmpdir  = "{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}",
        peerdone= "{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}/peer.done",
        origexpr= "{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}.txt.gz"
    output:
        outexpr = "{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}.txt.gz"
    run:
        import pandas as pd
        orig_df = pd.read_csv(f"{input[2]}", header=0, index_col=0, sep="\t")
        corr_df = pd.read_csv(f"{input[0]}/residuals.csv", header=None, index_col=None, sep=",")
        corr_df.index = list(orig_df.index)
        corr_df.columns = list(orig_df.columns)
        corr_df.to_csv(output[0], header=True, index=True, sep="\t")

rule calc_knn:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}.txt.gz"
    output:
        outexpr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_knn{knn}.txt.gz",
        outdist="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_knn{knn}.distmat.txt.gz"
    shell:
        "{config[python]} ../scripts/calc_knn.py --gx {input.expr} --k {wildcards.knn} --outdm {output.outdist} --outgx {output.outexpr}"