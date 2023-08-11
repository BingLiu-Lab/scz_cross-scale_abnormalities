library(tidyverse)
library(ggcorrplot)
library(Cairo)
library(fgsea)
library(gridExtra)
library(grid)
library(RRHO)
library(RColorBrewer)
library(readxl)
library(data.table)


base_dir = '/home/mwang/corticocortical/cell_transcriptional_enrichment_lake'
# base_dir = '/home/mwang/corticostriatal/cell_transcriptional_enrichment_lake'


# load dfc expr 
load(verbose=T, paste0(base_dir, '/lake_dfc_diff_expr.Rdata'))
# load vis expr
load(verbose=T, paste0(base_dir, '/lake_vis_diff_expr.Rdata'))

# combine dfc and vis expr
both_lake_cells = intersect(names(dfc_diff_expr_list), names(vis_diff_expr_list))
both_lake_sets  = list()
for (cell in both_lake_cells){
    # dfc
    dfc_dexpr_dat  = dfc_diff_expr_list[[cell]]
    dfc_dexpr_dat  = dfc_dexpr_dat[which(dfc_dexpr_dat$avg_logFC > 0),]
    dfc_dexpr_dat  = dfc_dexpr_dat[which(dfc_dexpr_dat$p_val_adj < .05),]

    # vis
    vis_dexpr_dat  = vis_diff_expr_list[[cell]]
    vis_dexpr_dat  = vis_dexpr_dat[which(vis_dexpr_dat$avg_logFC > 0),]
    vis_dexpr_dat  = vis_dexpr_dat[which(vis_dexpr_dat$p_val_adj < .05),]

    cell_genes = intersect(rownames(dfc_dexpr_dat), rownames(vis_dexpr_dat))
    print(cell)
    print(length(cell_genes))
    both_lake_sets[[cell]] = cell_genes
}

# gene markers for each cell
max_gene     = max(unlist(lapply(both_lake_sets, length)))
cell_gene_df = as.data.frame(matrix('', max_gene, length(both_lake_sets)))
colnames(cell_gene_df) = names(both_lake_sets)
for (cell in names(both_lake_sets)){
    write(cell,'')
    cell_genes  = both_lake_sets[[cell]]
    write_genes = c(cell_genes, rep('', max_gene - length(cell_genes)))
    cell_gene_df[cell] = write_genes
}
write_csv(cell_gene_df, paste0(base_dir, "/lake_cell_marker_genes.csv"))


# Run FGSEA for mean transcriptional correlations across datasets
gene_rank = read_excel(paste0(base_dir, '/corr_files/corticocortical_tmap_ahba_corr_avg.xlsx'))
# gene_rank = read_excel(paste0(base_dir, '/corr_files/corticostriatal_tmap_ahba_corr_avg.xlsx'))

gene_rank_vals = gene_rank$mean_r_vals
names(gene_rank_vals) = gene_rank$gene_name

fgseaRes = fgseaMultilevel(pathways = both_lake_sets,
                           stats = gene_rank_vals,
                           minSize = 15,
                           maxSize = 1000,
                           eps = 0,
                           nPermSimple = 10000)
fwrite(fgseaRes, file=paste0(base_dir, "/fgsea_lake_corticocortical/fgseaRes_mean.txt"), sep="\t", sep2=c("", " ", ""))
# fwrite(fgseaRes, file=paste0(base_dir, "/fgsea_lake_corticostriatal/fgseaRes_mean.txt"), sep="\t", sep2=c("", " ", ""))


# FGSEA enrichment score plots
# corticocortical
plotEnrichment(both_lake_sets[["Ex1"]], gene_rank_vals) + labs(title="corticocortical_lake_Ex1")
plotEnrichment(both_lake_sets[["Oli"]], gene_rank_vals) + labs(title="corticocortical_lake_Oli")

# corticostriatal
# plotEnrichment(both_lake_sets[["Ex4"]], gene_rank_vals) + labs(title="corticostriatal_lake_Ex4")
# plotEnrichment(both_lake_sets[["Oli"]], gene_rank_vals) + labs(title="corticostriatal_lake_Oli")