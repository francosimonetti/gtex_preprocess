#Define argv
#suppressMessages(library(tidyverse))

argv <- commandArgs(trailingOnly = TRUE)

# argv[1] <- '/data/malonso/GEU462_RPKM_QC_no_corrections_t.csv'
# argv[2] <- '/data/malonso/LDpruning/pca_geuv_allchr_3pc_t.csv'

exp_file = argv[1]
cov_file = argv[2]
out_file = argv[3]

message("Reading expression..")
gene_exp = read.table(file = exp_file, header = TRUE, sep = "\t", row.names = 1)


gene_exp_corr = gene_exp
#Read Covariates
pca = read.csv(file = cov_file, header = TRUE, sep="\t", row.names = 1)

if (nrow(pca) > 0) {
  message("Calculating cclm..")

  if (all(colnames(pca) %in% colnames(gene_exp)) & ncol(gene_exp) == ncol(pca))
  {n <- colnames(gene_exp)
  pca_ord <- pca[n]}

  if (! all(colnames(pca_ord) == colnames(gene_exp)))
  {stop("Rows and columns names don't match")}

  #Run our multiple linear regression and set the matrix of the residuals as our new expressions.
  #Making a copy of the gene expression transposed dataframe so that we can replace the values.

  #This loops through all the columns of the transposed gene expression which correspond
  #to each gene and for each gene it runs linear regression on the PEER factor covariates.
  #Then it sets the residuals to the new expression for that gene.
  for (i in 1:length(rownames(gene_exp))) {
    fit = lm(t(gene_exp[i,]) ~ t(as.matrix(pca_ord))) #CHEQUEAR dim peer fact y trasnsp
    expr_resid <- fit$residuals
    gene_exp_corr[i,] <- scale(expr_resid, center = TRUE, scale = TRUE)
  }
}

message("Writing results..")
#Write out the final expression file
column_names = colnames(gene_exp)
new_colnames = gsub("\\.", "-", column_names)
colnames(gene_exp_corr) <- new_colnames
if (grepl('.gz', out_file)) {
  tmpoutfile = gsub('.gz', '', out_file)
  write.table(gene_exp_corr , file = tmpoutfile, sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE)
  system(paste(sep=" ", 'gzip', tmpoutfile))
} else {
  write.table(gene_exp_corr , file = out_file, sep = "\t", row.names = TRUE, col.names=TRUE, quote=FALSE)
}
