#!/bin/bash

magma='/home/mwang/bin/magma'
base_dir=/home/mwang

eas_raw=${base_dir}/gene_analysis_out/magma_daner_natgen_pgc_eas.genes.raw
eur_raw=${base_dir}/gene_analysis_out/magma_daner_natgen_pgc_eur.genes.raw
meta_out=${base_dir}/gene_analysis_out/magma_daner_natgen_pgc_eas_eur_meta

# MAGMA meta gene results
${magma} --meta raw=${eas_raw}, ${eur_raw} --out ${meta_out}