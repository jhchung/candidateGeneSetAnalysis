options(stringsAsFactors = FALSE)
library(dplyr)
library(optparse)

# Define helper functions ----
gene_count_table <- function(data, gene_col = "Gene.refGene"){
  gene_table <- table(data[, gene_col]) %>%
    as.data.frame() %>%
    dplyr::rename(Gene  = Var1,
                  Count = Freq) %>%
    dplyr::arrange(desc(Count))
  return(gene_table)
}

# Define command line options ----
option_list <- list(
  make_option("--exonic",
              default = "C:/Lab/projects/WES/results/wgs100/chr22q11.2_region/22q_batch1.exonic.maf_all.withcounts",
              help = "Path to exonic variant file"),
  make_option("--deleterious",
              default = "C:/Lab/projects/WES/results/wgs100/chr22q11.2_region/22q_batch1.cadd_phred.maf_all.withcounts",
              help = "Path to deleterious variant file"),
  make_option("--deleterious_rare",
              default = "C:/Lab/projects/WES/results/wgs100/chr22q11.2_region/22q_batch1.cadd_phred.maf_0.01.withcounts",
              help = "Path to rare deleterious variant file"),
  make_option("--deleterious_novel",
              default = "C:/Lab/projects/WES/results/wgs100/chr22q11.2_region/22q_batch1.cadd_phred.maf_novel.withcounts",
              help = "Path to novel deleterious variant file"),
  make_option("--output_file",
              default = "C:/Lab/projects/WES/results/wgs100/chr22q11.2_region/22q_batch1.chr22q11.2_summary.txt",
              help = "Output directory for chr22q11.2 variant counts.")
)

args <- parse_args(OptionParser(option_list = option_list))

# Import data ----
exonic_var    <- read.delim(args$exonic, header = TRUE, sep = "\t")
del_var       <- read.delim(args$deleterious, header = TRUE, sep = "\t")
del_rare_var  <- read.delim(args$deleterious_rare, header = TRUE, sep = "\t")
del_novel_var <- read.delim(args$deleterious_novel, header = TRUE, sep = "\t")

# Count genes with variants ----
exonic_genes     <- gene_count_table(exonic_var)
del_genes        <- gene_count_table(del_var)
del_rare_genes   <- gene_count_table(del_rare_var)
del_novel_genes  <- gene_count_table(del_novel_var)

# Compile table ----
chr22_summary_table <- data.frame(
  "Variant category" = c("All exonic variants", "Predicted deleterious", 
                         "Rare deleterious", "Novel deleterious"),
  "Variants" = c(nrow(exonic_var), nrow(del_var), nrow(del_rare_var),
                 nrow(del_novel_var)),
  "Genes" = c(nrow(exonic_genes), nrow(del_genes), nrow(del_rare_genes),
              nrow(del_novel_genes))
)

# Export results ----
write.table(chr22_summary_table,
            args$output_file,
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
