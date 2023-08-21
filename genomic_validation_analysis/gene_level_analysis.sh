#!/bin/bash

magma='/home/mwang/bin/magma'
base_dir=/home/mwang

# EAS GWAS related files
gwas_snps_pval=${base_dir}/gwas_data_processed/gwas_snp_p/magma_daner_natgen_pgc_eas_snp_p.txt
gene_annot=${base_dir}/annot_out/magma_daner_natgen_pgc_eas.genes.annot
gene_analysis_out=${base_dir}/gene_analysis_out/magma_daner_natgen_pgc_eas
ref_bfile=${base_dir}/reference_data/g1000_eas/g1000_eas

# EUR GWAS related files
gwas_snps_pval=${base_dir}/gwas_data_processed/gwas_snp_p/magma_daner_natgen_pgc_eur_snp_p.txt
gene_annot=${base_dir}/annot_out/magma_daner_natgen_pgc_eur.genes.annot
gene_analysis_out=${base_dir}/gene_analysis_out/magma_daner_natgen_pgc_eur
ref_bfile=${base_dir}/reference_data/g1000_eur/g1000_eur

# MAGMA association test
${magma} --bfile ${ref_bfile} \
            --pval ${gwas_snps_pval}  ncol=N \
            --gene-annot ${gene_annot} \
            --out ${gene_analysis_out}