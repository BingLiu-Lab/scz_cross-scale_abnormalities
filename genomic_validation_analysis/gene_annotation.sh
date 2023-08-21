#!/bin/bash

magma='/home/mwang/bin/magma'
base_dir=/home/mwang
genome_build=${base_dir}/gene_location/NCBI37.3/NCBI37.3.gene.loc

# EAS GWAS related files
gwas_snps=${base_dir}/gwas_data_processed/gwas_snp_loc/magma_daner_natgen_pgc_eas_snp_loc.txt
annot_out=${base_dir}/annot_out/magma_daner_natgen_pgc_eas
ref_bfile=${base_dir}/reference_data/g1000_eas/g1000_eas

# EUR GWAS related files
gwas_snps=${base_dir}/gwas_data_processed/gwas_snp_loc/magma_daner_natgen_pgc_eur_snp_loc.txt
annot_out=${base_dir}/annot_out/magma_daner_natgen_pgc_eur
ref_bfile=${base_dir}/reference_data/g1000_eur/g1000_eur

# MAGMA gene annotation
${magma} --annotate window=5,5 \
            --snp-loc ${gwas_snps} \
            --gene-loc ${genome_build} \
            --out ${annot_out}