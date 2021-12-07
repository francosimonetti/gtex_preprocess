import sys
sys.path.append("../../knn_correction")
from utils.KNN import KNNcorrection as KNN
import numpy as np
import pandas as pd
import os
import argparse


def parse_args():

    parser = argparse.ArgumentParser("Script to calculate corrected KNN GX and distance matrix")

    parser.add_argument('--gx',
                        help="input gene expression file (.txt or .gz)",
                        type=str,
                        dest='gx_file')

    parser.add_argument('--dim',
                        help="number of PCA dimensions for sample distance (between 0 and 1 as fraction of nsamples, or > 1 for a specific dimension",
                        type=float,
                        default=None,
                        dest='pcadim')

    parser.add_argument('--k',
                        help="k nearest neighbors",
                        type=int,
                        dest='k')

    parser.add_argument('--outdm',
                        help="output distance matrix",
                        type=str,
                        default=None,
                        dest='outdm')

    parser.add_argument('--outgx',
                        help="output file for corrected gx",
                        type=str,
                        dest='outgx')

    opts = parser.parse_args()
    return opts

if __name__ == "__main__":
    
    args = parse_args()

    expr_df = pd.read_csv(args.gx_file, header=0, index_col=0, sep="\t")
    samplenames = list(expr_df.columns)
    knn = KNN(expr_df.values.T, dim = args.pcadim)
    gx_corr, gx_term = knn.knn_corr(K = args.k)

    gxcorr_df = pd.DataFrame(gx_corr.T, columns=samplenames, index=expr_df.index)
    gxcorr_df.to_csv(args.outgx, sep="\t", header=True, index=True)

    if args.outdm is not None:
        distmat = pd.DataFrame(knn.dm, columns=samplenames, index=samplenames)
        distmat.to_csv(args.outdm, sep="\t", header=True, index=True)

