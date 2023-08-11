library(tidyverse)
library(ggcorrplot)
library(Cairo)
library(fgsea)
library(gridExtra)
library(grid)
library(readxl)
library(data.table)


base_dir = '/home/mwang/corticocortical/cell_transcriptional_enrichment_aba'
# base_dir = '/home/mwang/corticostriatal/cell_transcriptional_enrichment_aba'


# load mtg expr
load(verbose=T, paste0(base_dir, '/aba_mtg_diff_expr.Rdata'))

fgsea_mtg = NULL
for (cell in names(mtg_diff_expr_list)){
  cell_dexpr     = mtg_diff_expr_list[[cell]]
  cell_dexpr_df  = cell_dexpr[which(cell_dexpr$avg_logFC > 0),]
  cell_dexpr_df  = cell_dexpr_df[which(cell_dexpr_df$p_val_adj < .05),]
  
  fgsea_mtg[[cell]] = rownames(cell_dexpr_df)
}

# gene markers for each cell
max_gene     = max(unlist(lapply(fgsea_mtg, length)))
cell_gene_df = as.data.frame(matrix('', max_gene, length(fgsea_mtg)))
colnames(cell_gene_df) = names(fgsea_mtg)
for (cell in names(fgsea_mtg)){
    write(cell,'')
    cell_genes  = fgsea_mtg[[cell]]
    write_genes = c(cell_genes, rep('', max_gene - length(cell_genes)))
    cell_gene_df[cell] = write_genes
}
write_csv(cell_gene_df, paste0(base_dir, "/mtg_cell_marker_genes.csv"))


# Run FGSEA for mean transcriptional correlations across datasets
gene_rank = read_excel(paste0(base_dir, '/corr_files/corticocortical_tmap_ahba_corr_avg.xlsx'))
# gene_rank = read_excel(paste0(base_dir, '/corr_files/corticostriatal_tmap_ahba_corr_avg.xlsx'))

gene_rank_vals = gene_rank$mean_r_vals
names(gene_rank_vals) = gene_rank$gene_name

fgseaRes = fgseaMultilevel(pathways = fgsea_mtg,
                           stats = gene_rank_vals,
                           minSize = 15,
                           maxSize = 1000,
                           eps = 0,
                           nPermSimple = 10000)
fwrite(fgseaRes, file=paste0(base_dir, "/fgsea_aba_corticocortical/fgseaRes_mean.txt"), sep="\t", sep2=c("", " ", ""))
# fwrite(fgseaRes, file=paste0(base_dir, "/fgsea_aba_corticostriatal/fgseaRes_mean.txt"), sep="\t", sep2=c("", " ", ""))


# FGSEA enrichment score plots
# corticocortical
plotEnrichment(fgsea_mtg[["Exc_Ex1(ABA)"]], gene_rank_vals) + labs(title="corticocortical_aba_Exc_Ex1(ABA)")
plotEnrichment(fgsea_mtg[["Oligo_OPALIN"]], gene_rank_vals) + labs(title="corticocortical_aba_Oligo_OPALIN")

# corticostriatal
# plotEnrichment(fgsea_mtg[["Exc_RORB"]], gene_rank_vals) + labs(title="corticocortical_aba_Exc_RORB")
# plotEnrichment(fgsea_mtg[["Oligo_OPALIN"]], gene_rank_vals) + labs(title="corticocortical_aba_Oligo_OPALIN")