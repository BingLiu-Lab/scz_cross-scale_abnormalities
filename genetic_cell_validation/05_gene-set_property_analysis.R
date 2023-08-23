library(EWCE)
library(MAGMA.Celltyping)
library(Seurat)
library(readr)
data(all_hgnc_wtEntrez)


# Cell specific gene expression from Lake dfc (NOTE: the analysis is similar for Lake vis and ABA mtg)
base_dir = '/home/mwang/cell_type_data/lake_dfc'
load(paste0(base_dir, '/cell_type_ewce_lake_dfc.rda'), verbose=T)
dfc_ctd = ctd

fine_cells = colnames(dfc_ctd[[1]]$mean_exp)
dfc_fine_mean_exp         = as.data.frame(dfc_ctd[[1]]$mean_exp)
dfc_fine_mean_exp$Average = rowMeans(dfc_fine_mean_exp)
dfc_fine_mean_exp$hgnc    = rownames(dfc_fine_mean_exp)

# map hgnc to entrez id and format data for write
dfc_fine_mean_exp_name = merge(x=dfc_fine_mean_exp, y=all_hgnc_wtEntrez, by.x='hgnc', by.y='hgnc_symbol')
dfc_fine_gene_out   = dfc_fine_mean_exp_name[c('entrezgene', fine_cells, 'Average')]

save_dir = '/home/mwang/gene_property_analysis_out_lake_dfc'
gene_cond_file = paste0(save_dir, '/lake_dfc_expr_property.txt')
write_delim(x=dfc_fine_gene_out, gene_cond_file, delim='\t')


# run gene property analysis in terminal
magma --gene-results /home/mwang/gene_analysis_out/magma_daner_natgen_pgc_eas_eur_meta.genes.raw \
--gene-covar /home/mwang/gene_property_analysis_out_lake_dfc/lake_dfc_expr_property.txt \
--model direction=pos condition='Average' \
--out /home/mwang/gene_property_analysis_out_lake_dfc/lake_dfc_expr_property_result