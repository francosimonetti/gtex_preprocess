import sys
sys.path.append("/data/franco/knn_correction/utils")
from KNN import KNNcorrection as KNN
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

    parser.add_argument('--gtpca',
                        help="input genotype PCs (I x P)",
                        type=str,
                        default=None,
                        dest='gtpca_file')

    parser.add_argument('--dim',
                        help="number of PCA dimensions for sample distance (between 0 and 1 as fraction of nsamples, or > 1 for a specific dimension",
                        type=float,
                        default=None,
                        dest='pcadim')

    parser.add_argument('--k',
                        help="k nearest neighbors",
                        type=int,
                        dest='k')

    parser.add_argument('--kreciprocal',
                        help="use k-reciprocal nearest neighbors",
                        action='store_true',
                        default=False,
                        dest='use_kreciprocal')

    parser.add_argument('--outdm',
                        help="output distance matrix",
                        type=str,
                        default=None,
                        dest='outdm')

    parser.add_argument('--outgx',
                        help="output file for corrected gx",
                        type=str,
                        dest='outgx')

    # parser.add_argument('--is-transposed',
    #                     help="transpose input expr matrix (GxN)",
    #                     action='store_true',
    #                     default=False,
    #                     dest='is_transposed')

    opts = parser.parse_args()
    return opts

if __name__ == "__main__":
    
    args = parse_args()

    filedir = os.path.dirname(args.outgx)
    if not os.path.exists(filedir) and filedir != "":
        os.makedirs(filedir)

    if args.outdm is not None:
        filedir = os.path.dirname(args.outdm)
        if not os.path.exists(filedir) and filedir != "":
            os.makedirs(filedir)

    expr_df = pd.read_csv(args.gx_file, header=0, index_col=0, sep="\t")
    genenames = list(expr_df.columns)
    samplenames = list(expr_df.index)

    if args.gtpca_file is not None:
        gtpca_df = pd.read_csv(args.gtpca_file, header=0, index_col=0, sep="\t")
        if gtpca_df.shape[0] < gtpca_df.shape[1]:
            print("More PCs than samples? Transposing matrix")
            gtpca_df = gtpca_df.T
        gt_samplenames = list(gtpca_df.index)
        common_samples = [x for x in gt_samplenames if x in samplenames]
        if len(common_samples) != len(samplenames):
            raise Exception("Not all samples have genotype PCs. Please check that samples match")
        gtpca_sorted_df = gtpca_df.loc[samplenames]

    if args.k > 0:
        knn = KNN(expr_df.values, dim = args.pcadim)
        if args.gtpca_file is not None:
            gx_corr, gx_term = knn.knn_corr(K = args.k, dmdata=gtpca_sorted_df.values, dim = args.pcadim, kr=args.use_kreciprocal)
        else:
            gx_corr, gx_term = knn.knn_corr(K = args.k, kr=args.use_kreciprocal, dim = args.pcadim)
        gxcorr_df = pd.DataFrame(gx_corr, columns=genenames, index=expr_df.index)
        gxcorr_df.to_csv(args.outgx, sep="\t", header=True, index=True)
        if args.outdm is not None:
            distmat = pd.DataFrame(knn.dm, columns=samplenames, index=samplenames)
            distmat.to_csv(args.outdm, sep="\t", header=True, index=True)
    else:
        expr_df.to_csv(args.outgx, sep="\t", header=True, index=True)
