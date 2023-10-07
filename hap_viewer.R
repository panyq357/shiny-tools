
config = list(
    "vcf" = "/mnt/d/Analysis/YXH.2023-04-24.GWAS/resources/ALL_GW.vcf.gz"
)

get_df_from_vcf = function(vcf_path, chr, start, end) {
    # Given a region, subset vcf, and read into a data.frame.
    vcf_content <- system2(
        command = "tabix",
        args = c("-h", vcf_path, paste0(chr, ":", start, "-", end)),
        stdout = T
    )
    vcf_content = vcf_content[grep("^##", vcf_content, invert=T)]
    df = read.table(text=paste0(vcf_content, collapse="\n"), header=T, comment.char="", check.names=F)
    return(df)
}

var_stats = function(vcf_df) {

    row.names(vcf_df) = vcf_df$ID

    out = list()
    out$info = vcf_df[c("#CHROM", "POS", "REF", "ALT")]
    out$genotype = as.matrix(vcf_df[10:length(vcf_df)])

    for (i in 1:nrow(out$info)) {
        
        ref = out$info[i, "REF"]
        alt = out$info[i, "ALT"]

        out$info[i, "REF_FREQ"] = mean(out$genotype[i,] == "0/0")
        out$info[i, "ALT_FREQ"] = mean(out$genotype[i,] == "1/1")

        out$genotype[i,][(out$genotype[i,] != "0/0") & (out$genotype[i,] != "1/1")] = "N"
        out$genotype[i,][out$genotype[i,] == "0/0"] = ref
        out$genotype[i,][out$genotype[i,] == "1/1"] = alt
    }

    return(out)
}

# out = var_stats(vcf_df)

divide_haps = function(genotype, var_id=NULL, min_n=10, min_maf=0.2, min_missing=0.1) {
    # Out: a list, each element is a vector containing sample ID of that haplotype.

    hap_list = list()

    for (i in 1:ncol(genotype)) {
        hap = paste0(genotype[,i], collapse="")
        id = colnames(genotype)[i]

        if (! length(hap_list[[hap]]) > 0) {
            hap_list[[hap]] = list(id)
        } else {
            hap_list[[hap]] = append(hap_list[[hap]], id)
        }
    }

    return(hap_list)
}

# chr = 3
# start = 16729501
# end = 16735109
# vcf_df = get_df_from_vcf(config[["vcf"]], chr, start, end)

