

import json
with open("gtex_v8.custom.config") as instr:
    gtex_conf = json.load(instr)

configfile: "config.json",

target_tissues = ['ebv'] #['as', 'ebv', 'ms', 'wb'] #gtex_conf['tshorts']
target_gxnorm  = ["tpms_qcfilter"]

# subworkflow part2:
#     snakefile:
#         "pipe2.snakefile"

PCA   = [0]
NPEER = [0]
K     = [10, 20, 30]
DIM   = ["1"]
KNN   = ["knn"]
GTDIM = ["3"]

# GTKNN = ["knn"]
# GTK   = [0]
# GTDIM = ["0"]

wildcard_constraints:
    tissue  = "|".join([x for x in target_tissues]),
    gxnorm = "|".join([x for x in target_gxnorm]),
    npeer  = "\d+",
    pca    = "\d+",
    k      = "\d+",
    knn    = "|".join([x for x in KNN]),
    gtdim = "|".join([x for x in GTDIM]),

### To obtain peer numbers according to GTEx
# tissue_peer_dict = dict()
# with open(config['tissue_table']) as ifile:
#     for line in ifile:
#         if line.startswith("#"): continue
#         arr = line.rstrip().split("\t")
#         tissue_peer_dict[arr[1]] = arr[3]

rule all:
    input:
        #expand(config['covariate_dir']+"/{tissue}.pca.txt", tissue=target_tissues, pca=PCA)
        #expand("{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}.txt.gz", tissue=target_tissues, pca=PCA, gxnorm=target_gxnorm, workdir=config['gtex_exprdir'])
        #expand("{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}.txt.gz", tissue=target_tissues, pca=PCA, gxnorm=target_gxnorm, npeer=NPEER, workdir=config['gtex_exprdir'])
        #expand("{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_{knn}{k}d{dim}.txt.gz", tissue=target_tissues, pca=PCA, gxnorm=target_gxnorm, npeer=NPEER, k=K, knn=KNN, dim=DIM, workdir=config['gtex_exprdir'])
        expand("{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_{knn}GT{k}d{gtdim}_{knn}{k}d{dim}.txt.gz", tissue=target_tissues, pca=PCA, gxnorm=target_gxnorm, npeer=NPEER, k=K, knn=KNN, dim=DIM, gtdim=GTDIM, workdir=config['gtex_exprdir'])
        

def get_covariate_input(wildcards):
    return f"{config['covariate_dir']}/{gtex_conf[wildcards.tissue]['fullname']}.v8.covariates.txt"

rule prepare_pca_covariates:
    input:
        get_covariate_input
    output:
        outfile=config['covariate_dir']+"/{tissue}.pca.txt",
    shell:
        #"grep -iv inferred {input} > {output.outfile}.tmp ; head -n $(( {wildcards.pca} + 1 )) {output.outfile}.tmp > {output.outfile}; rm {output.outfile}.tmp"
        "head -n 1 {input} > {output.outfile}; sleep 6;"
        "tail -n +2 {input} | grep 'PC' >> {output.outfile}"
       
rule cclm_correction:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.txt.gz",
        cov=config['covariate_dir']+"/{tissue}.pca.txt"
    output:
        "{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}.txt.gz"
    shell:
        "Rscript ../scripts/cclm.R {input.expr} {input.cov} {output} {wildcards.pca}"

# decompress expr, sort covariates as expr samples, transpose covariate, delete row and col names for covariate
rule prepare_peer:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}.txt.gz",
        # cov=config['covariate_dir']+"/{tissue}.pca{pca}.txt"
    output:
        tmpdir=directory("{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}")
    shell:
        "{config[python]} ../scripts/prepare_peer_inputs.py --gx {input.expr} --outdir {output.tmpdir}" #--cov {input.cov} 

## Peer with covariate file cov.tab
# rule calc_peer:
#     input:
#         tmpdir = "{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}"
#     output:
#         touch("{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}/peer.done")
#     shell:
#         "if [ {wildcards.npeer} -ne 0 ] ; then peertool -f {input.tmpdir}/gx.tab -n {wildcards.npeer} -c {input.tmpdir}/cov.tab --no_a_out --no_z_out --e_pb 10 --e_pa 0.1 --a_pb 0.01 --a_pa 0.001 -o {input.tmpdir};"
#         "else sed 's/\t/,/g' {input.tmpdir}/gx.tab > {input.tmpdir}/residuals_t.csv; awk -f ../scripts/transpose_csv.awk {input.tmpdir}/residuals_t.csv >  {input.tmpdir}/residuals.csv; fi;"

rule calc_peer_nocov:
    input:
        tmpdir = "{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}"
    output:
        touch("{workdir}/{gxnorm}/tmp_{tissue}_pca{pca}_npeer{npeer}/peer.done")
    shell:
        ### PEER with input covariates or with previous cclm of those covariates is the same! yay!
        "if [ {wildcards.npeer} -ne 0 ] ; then peertool -f {input.tmpdir}/gx.tab -n {wildcards.npeer} --no_a_out --no_z_out --e_pb 10 --e_pa 0.1 --a_pb 0.01 --a_pa 0.001 -o {input.tmpdir};"
        "else sed 's/\t/,/g' {input.tmpdir}/gx.tab > {input.tmpdir}/residuals_t.csv; awk -f ../scripts/transpose_csv.awk {input.tmpdir}/residuals_t.csv >  {input.tmpdir}/residuals.csv; fi;"

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
        corr_df.T.to_csv(output[0], header=True, index=True, sep="\t") # transposed matrix here!

rule calc_knn:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}.txt.gz"
    output:
        outexpr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_{knn}{k}d{dim}.txt.gz",
    shell:
        "if [ '{wildcards.knn}' = 'knn' ] ; then {config[python]} ../scripts/calc_krnn.py --gx {input.expr} --k {wildcards.k} --dim {wildcards.dim} --outgx {output.outexpr}; elif [ '{wildcards.knn}' = 'krnn' ] ; then {config[python]} ../scripts/calc_krnn.py --gx {input.expr} --k {wildcards.k} --dim {wildcards.dim} --outgx {output.outexpr} --kreciprocal; fi"

rule calc_GTknn:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}.txt.gz",
        gtpca=config['covariate_dir']+"/{tissue}.pca.txt"
    output:
        outexpr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_{knn}GT{k}d{gtdim}.txt.gz"
    shell:
        "if [ '{wildcards.knn}' = 'knn' ] ; then {config[python]} ../scripts/calc_krnn.py --gx {input.expr} --gtpca {input.gtpca} --k {wildcards.k} --dim {wildcards.gtdim} --outgx {output.outexpr}; elif [ '{wildcards.knn}' = 'krnn' ] ; then {config[python]} ../scripts/calc_krnn.py --gx {input.expr} --gtpca {input.gtpca} --k {wildcards.k} --dim {wildcards.gtdim} --outgx {output.outexpr} --kreciprocal; fi"

rule calc_double_knn:
    input:
        expr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_{knn}GT{k}d{gtdim}.txt.gz",
    output:
        outexpr="{workdir}/{gxnorm}/{tissue}.{gxnorm}.pc_lncrna.pca{pca}_npeer{npeer}_{knn}GT{k}d{gtdim}_{knn}{K}d{dim}.txt.gz"
    shell:
        "if [ '{wildcards.knn}' = 'knn' ] ; then {config[python]} ../scripts/calc_krnn.py --gx {input.expr} --k {wildcards.k} --dim {wildcards.dim} --outgx {output.outexpr}; elif [ '{wildcards.knn}' = 'krnn' ] ; then {config[python]} ../scripts/calc_krnn.py --gx {input.expr} --k {wildcards.k} --dim {wildcards.dim} --outgx {output.outexpr} --kreciprocal; fi"
