library(GenomicFeatures)

# Prepare txdb for gene annotation track.
gene_gr <- rtracklayer::import("/home/panyq/Tools/index-scripts/os/rap-db/2023-10-25/gtf/os.rap-db.make.gtf")
seqlevels(gene_gr) <- sub("chr[0]*", "", seqlevels(gene_gr))
txdb <- makeTxDbFromGRanges(gene_gr)

# Prepare gene annotation data.frame.
gene_anno <- readr::read_tsv("/home/panyq/Tools/index-scripts/os/rap-db/2023-10-25/rawdata/IRGSP-1.0_representative_annotation_2023-09-07.tsv.gz")
gene_anno <- gene_anno[c("Locus_ID", "Description", "Oryzabase Gene Symbol Synonym(s)")]
gene_anno <- gene_anno[!duplicated(gene_anno$Locus_ID),]

# Prepare GWAS result data.frame list.
result_list <- list(
    "Blink_GL" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/GL/GAPIT.Blink.GL.GWAS.Results.csv",
    "CMLM_GL" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/GL/GAPIT.CMLM.GL.GWAS.Results.csv",
    "GLM_GL" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/GL/GAPIT.GLM.GL.GWAS.Results.csv",
    "Blink_GW" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/GW/GAPIT.Blink.GW.GWAS.Results.csv",
    "CMLM_GW" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/GW/GAPIT.CMLM.GW.GWAS.Results.csv",
    "GLM_GW" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/GW/GAPIT.GLM.GW.GWAS.Results.csv",
    "Blink_TGW" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/TGW/GAPIT.Blink.TGW.GWAS.Results.csv",
    "CMLM_TGW" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/TGW/GAPIT.CMLM.TGW.GWAS.Results.csv",
    "GLM_TGW" = "/mnt/d/Analysis/LYH.2022-04-01.3K-RG-2K-samples-GWAS/results/TGW/GAPIT.GLM.TGW.GWAS.Results.csv"
) |> lapply(function(x) {
    df <- readr::read_csv(x)
    df <- df[c("Chromosome", "Position", "FDR_Adjusted_P-values")]
    names(df) <- c("chr", "pos", "p")
    return(df)
})


config = list(
    title = "3K GWAS Viewer",

    # A GRanges for Gviz to plot gene structure.
    txdb = txdb,

    # A data.frame, first column is gene_id, all other columns will be shown in shiny app.
    gene_anno = gene_anno,

    # A data.frame list, items are data.frames with three columns: "chr", "pos" and "p".
    result_list = result_list
)

