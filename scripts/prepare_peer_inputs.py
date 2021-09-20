import numpy as np
import pandas as pd
import os
import argparse


def parse_args():

    parser = argparse.ArgumentParser("Script to prepare input files format for PEER")

    parser.add_argument('--gx',
                        help="input gene expression file compressed .gz",
                        type=str,
                        dest='gx_file')

    parser.add_argument('--cov',
                        help="input Covariate file",
                        type=str,
                        dest='cov_file')

    parser.add_argument('--outdir',
                        help="output file for gx",
                        type=str,
                        dest='outdir')

    # parser.add_argument('--outgx',
    #                     help="output file for gx",
    #                     type=str,
    #                     dest='outgx')

    # parser.add_argument('--outcov',
    #                     help="output file for covs",
    #                     type=str,
    #                     dest='outcov')

    opts = parser.parse_args()
    return opts

if __name__ == "__main__":
    
    args = parse_args()

    gx_df  = pd.read_csv(args.gx_file, header=0, index_col=0, sep="\t")
    cov_df = pd.read_csv(args.cov_file, header=0, index_col=0, sep="\t")

    if gx_df.shape[0] > gx_df.shape[1]: # it's GxN
        # transpose to NxG
        gx_df = gx_df.T
        gx_df.index.name = "gene_id"

    nsamples = gx_df.shape[0]
    if cov_df.shape[0] != cov_df.shape[1]:
        if cov_df.shape[0] != nsamples:
            cov_df = cov_df.T
    else:
        print("ERROR! Same number of covariates as samples!")
        raise

    if cov_df.shape[0] != nsamples:
        print("FATAL ERROR! covariates have different nsamples as gx")
        raise

    sorted_cov_df = cov_df.loc[list(gx_df.index)]

    # outdir = os.path.dirname(args.outgx)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    with open(os.path.join(args.outdir, "gx.tab"), 'w') as ofile:
        gx_df.to_csv(ofile, header=False, index=False, sep="\t")
    
    with open(os.path.join(args.outdir, "cov.txt"), 'w') as ofile:
        sorted_cov_df.to_csv(ofile, header=False, index=False, sep="\t")

    # with open(args.outgx, 'w') as ofile:
    #     gx_df.to_csv(ofile, header=False, index=False, sep="\t")
    
    # with open(args.outcov, 'w') as ofile:
    #     sorted_cov_df.to_csv(ofile, header=False, index=False, sep="\t")
    