library(readr)

# for EAS population
sz_gwas_path = "/home/mwang/gwas_data_from_PGC/daner_natgen_pgc_eas"
sz_gwas = read_delim(sz_gwas_path, delim=' ')
sz_gwas$N = sz_gwas$Nca + sz_gwas$Nco
 
sz_write = sz_gwas[c('SNP','CHR','BP')]
base_path = "/home/mwang/gwas_data_processed/gwas_snp_loc"
sz_write_path = paste0(base_path, '/magma_daner_natgen_pgc_eas_snp_loc.txt')
print(sz_write_path)
write_delim(sz_write, sz_write_path, delim='\t', col_names=FALSE)
 
sz_write = sz_gwas[c('SNP','CHR','BP','P','N')]
base_path = "/home/mwang/gwas_data_processed/gwas_snp_p"
sz_write_path = paste0(base_path, '/magma_daner_natgen_pgc_eas_snp_p.txt')
print(sz_write_path)
write_delim(sz_write, sz_write_path, delim='\t', col_names=TRUE)


# for EUR population
sz_gwas_path = "/home/mwang/gwas_data_from_PGC/daner_natgen_pgc_eur"
sz_gwas = read_delim(sz_gwas_path, delim=' ')
sz_gwas$N = sz_gwas$Nca + sz_gwas$Nco

sz_write = sz_gwas[c('SNP','CHR','BP')]
base_path = "/home/mwang/gwas_data_processed/gwas_snp_loc"
sz_write_path = paste0(base_path, '/magma_daner_natgen_pgc_eur_snp_loc.txt')
print(sz_write_path)
write_delim(sz_write, sz_write_path, delim='\t', col_names=FALSE)

sz_write = sz_gwas[c('SNP','CHR','BP','P','N')]
base_path = "/home/mwang/gwas_data_processed/gwas_snp_p"
sz_write_path = paste0(base_path, '/magma_daner_natgen_pgc_eur_snp_p.txt')
print(sz_write_path)
write_delim(sz_write, sz_write_path, delim='\t', col_names=TRUE)