import sys, os
#sys.path.append("..")
import readgtf
from collections import defaultdict
import pandas as pd
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Filters gene expression by genes and donors')

    parser.add_argument('--gx',
                         help='input file with gene expression (genes x donors)',
                         type=str,
                         dest='gx_file')

    parser.add_argument('--donors',
                         dest='donor_file',
                         type=str,
                         default=None,
                         help='file with donor ids')

    parser.add_argument('--out',
                        dest='out_file',
                        default=None,
                        help='output file for new gene expression')    

    parser.add_argument('--dataset',
                        dest='dataset',
                        default='gtex',
                        help='gtex or cardiogenics')

    parser.add_argument('--gtf',
                         type=str,
                         help='GTF file to use',
                         dest='gtf_file')

    parser.add_argument('--biotype',
                        nargs = '*',
                        dest = 'biotype',
                        type = str,
                        help = 'type of genes to select')

    opts = parser.parse_args()
    return opts

def read_samples(donorfile):    
    with open(donorfile, 'r') as samfile:
        sample = 0
        samplenames = list()
        # skip first two lines
        next(samfile)
        next(samfile)
        for line in samfile:
            if re.search('^#', line):
                continue
            sample += 1
            samplenames.append(line.strip().split()[0])
    return samplenames

def filter_donors(df, donors):
    if donors is not None:
        donor_list = df.columns
        common  = [x for x in donors if x in donor_list]
        print("{:d} donors remained from {:d}".format(len(common), len(donor_list)))
        return df[common]
    else:
        return df

def filter_rows(df, genedict):
    gx_gene_list = df.index
    common  = [genedict[x] for x in gx_gene_list]
    print("{:d} genes remained from {:d}".format(sum(common), len(gx_gene_list)))
    return df[common]

if __name__ == '__main__':
    args = parse_args()

    print("Gene Expr File: {:s}".format(args.gx_file))

    # can't read this with current library
    # gtfpath = "/cbscratch/franco/datasets/gtex/gencode.v28lift37.annotation.gtf.gz"
    gtfpath = args.gtf_file
    print("Reading GENCODE file")

    if args.donor_file is None:
        donors = None
    else:
        donors = read_samples(args.donor_file)

    if args.dataset == "gtex":
        gene_info = readgtf.gencode(gtfpath, trim=False, biotype=args.biotype)
    if args.dataset == "cardiogenics":
        gene_info = readgtf.gencode(gtfpath, trim=True, biotype=args.biotype)

    gene_dict = defaultdict(lambda: False)
    for g in gene_info:
        gene_dict[g.ensembl_id] = True 	

    gx_df = pd.read_table(args.gx_file, sep="\t", header=0, index_col=0)
    new_gx_df = filter_rows(gx_df, gene_dict)
    sorted_gx_df = filter_donors(new_gx_df, donors)
    # outfile = GXFILE+".gencode_filter"
    if args.out_file is None:
        args.out_file = args.gx_file + ".{:s}_filtered".format("_".join(args.biotype))
        print("Outfile:", args.out_file)
    outdir = os.path.dirname(os.path.realpath(args.out_file))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    sorted_gx_df.to_csv(args.out_file, doublequote=False, sep="\t")
