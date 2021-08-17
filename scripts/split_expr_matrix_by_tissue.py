#!/usr/bin/env python3
# Author: Saikat Banerjee and Franco Simonetti

import pandas as pd
import argparse
import gzip
import os

def parse_args():

    parser = argparse.ArgumentParser(description='Extract samples of a particular tissue from GTEx gene expression file')

    parser.add_argument('--input',
                        type=str,
                        dest='infilepath',
                        metavar='FILE',
                        help='input GCT file containing all samples')

    parser.add_argument('--outdir',
                        type=str,
                        dest='outdir',
                        metavar='FILE',
                        help='output directory for tissue expression files')

    parser.add_argument('--tissuefile',
                        type=str,
                        dest='tissue_table',
                        metavar='STR',
                        help='tissue table with matching acronyms and abbreviations') #"/data/franco/datasets/gtex_v8/tissues_list.txt"

    parser.add_argument('--pheno',
                        type=str,
                        dest='phenofilepath',
                        metavar='STR',
                        help='phenotype file of the GTEx consortium')

    parser.add_argument('--donors',
                        type=str,
                        dest='donor_file',
                        metavar='STR',
                        help='File with donor ids with genotype')

    opts = parser.parse_args()
    return opts


def get_samples(pheno_df, tissue, gt_donors):
    sub_df = pheno_df.loc[(pheno_df['SMTORMVE'] != "FLAGGED") & (pheno_df['SMGEBTCHT'] == "TruSeq.v1") & (pheno_df['SMTSD'] == tissue)]
    valid_donors = [i for i in list(sub_df["SAMPID"]) if '-'.join(i.split('-')[:2]) in gt_donors]
    # print([True for i in list(sub_df["SAMPID"]) if i in valid_donors])
    return valid_donors
    
def get_donors(path):
    donor_ids = list()
    with open(path, 'r') as instream:
        next(instream)
        next(instream)
        for line in instream:
            donor_ids.append(line.strip().split()[0])
    return donor_ids

def read_gct(gct_file, donor_ids=None):
    """
    Load GCT as DataFrame
    """    
    df = pd.read_csv(gct_file, sep='\t', skiprows=2, index_col=0)
    df.drop('Description', axis=1, inplace=True)
    df.index.name = 'gene_id'
    if donor_ids is not None:
        df = df[[i for i in df.columns if '-'.join(i.split('-')[:2]) in donor_ids]]
    return df

def write_gct(df, filepath, trim_ids=True, header=False):
    """
    Write dataframe as a GCT file
    """
    with open(filepath, 'w') as mfile:
        if header:
            mfile.write("#1.2\n")
            mfile.write('%i\t%i\n' % (df.shape[0], df.shape[1] - 1))
        # mfile.write(str(df.shape[0])+'\t'+str(df.shape[1] - 1)+'\n')
        if trim_ids:
            new_headers = ['-'.join(i.split('-')[:2]) for i in df.columns]
            df.columns = new_headers
        df.to_csv(mfile, sep='\t', index=True, header=True)
    print ('New GCT file written with %i samples and %i genes.\n' % (df.shape[1], df.shape[0] - 1))
    

if __name__=='__main__':

    opts = parse_args()

    # Load genotype donor ids, we will use them to filter rna-seq samples that also have genotype
    pheno_df = pd.read_csv(opts.phenofilepath, sep="\t", comment="#", header=0)
    gt_donors = get_donors(opts.donor_file)

    # Get list of all tissues present in the phenotype file
    tissue_names = list(pheno_df[pheno_df['SMGEBTCHT'].str.contains("TruSeq.v1")]["SMTSD"].unique())
    
    # create a dictionary from the whole name to the acronyms
    names_df = pd.read_table(opts.tissue_table, sep="\t", comment="#", header=None)
    names_df.columns = ["long_name", "acr", "id", "npeer"]
    names2acr = dict(zip(names_df["long_name"], names_df["acr"]))

    # From the phenotype file loaded above, filter RNA-seq samples that are now flagged and belong to each tissue    
    donors_per_tissue = dict()
    for tissue in tissue_names:
        donors_per_tissue[tissue] = get_samples(pheno_df, tissue, gt_donors)

    # Load expression or counts
    outtype = "tpm"
    if "reads" in opts.infilepath:
        outtype="reads"
    expression_df = read_gct(opts.infilepath, gt_donors)
    print ('Number of samples read from GCT file: %i\n' % (expression_df.shape[1] - 1))

    if not os.path.exists(opts.outdir):
        os.makedirs(opts.outdir)
    # Count for each tissue, the number of samples with Genotype, should match gtexportal
    # Filter tissues with less than 30 samples
    counter = 1
    for t in donors_per_tissue:
        if len(donors_per_tissue[t]) > 30 and t in names2acr:
            print(counter, t, len(donors_per_tissue[t]))
            counter += 1

            print(f"Processing {t} -> {names2acr[t]}")
            outfile = os.path.join(opts.outdir, f"GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9.gene_{outtype}.{names2acr[t]}.txt")
            if not os.path.exists(outfile):
                expr_tissue_df = expression_df[donors_per_tissue[t]]
                write_gct(expr_tissue_df, outfile)
            else:
                print(f"{t} file exists")
            # compress the file
            with open(outfile, 'rb') as f_in, gzip.open(outfile+'.gz', 'wb') as f_out:
                f_out.writelines(f_in)
        else:
            print(f"ERROR: tissue {t} not in table or <30 samples")