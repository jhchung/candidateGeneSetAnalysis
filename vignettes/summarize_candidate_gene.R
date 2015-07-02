#' Identify variants in genes found in the target gene set and Count 
#' case/control numbers .
#' Creates Figures 

# Setup ----
options(stringsAsFactors = FALSE)
library(tbx1GeneNetwork) # Personal package
library(optparse)
library(plyr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(foreach)
library(doSNOW)

################################################################################
# Begin script 
################################################################################
# Define command line options ----
option_list <- list(
  make_option("--cohort1",
              default =
                "C:/Lab/projects/WES/results/wes186/variant_filtering/22q_exomes.cadd_phred.maf_0.01.withcounts",
              help = "Path to exonic variant file"),
  make_option("--cohort2",
              default ="C:/Lab/projects/WES/results/wgs100/variant_filtering/22q_batch1.cadd_phred.maf_0.01.withcounts",
              help = "Path to variant file for cohort 2"),
  make_option("--gmt",
              default = "C:/Lab/projects/WES/data/gene_lists/pax9_stringdb.gmt",
              help = "Path to gene set file."),
  make_option("--gene_set",
              default = "PAX9_interaction_network",
              help = "Name of target gene set"),
  make_option("--cohort1_pheno",
              default = "C:/Lab/projects/WES/data/wes_samples_phenotypes/wes_updated_phenotype_4_11_15.txt",
              help = "Phenotype data for cohort 1 in plink format"),
  make_option("--cohort2_pheno",
              default = "C:/Lab/projects/WES/data/wgs100/wgs100_phenotype-01_16_15.txt",
              help = "Phenotype data for cohort 2 in plink format"),
  make_option("--cohort1_alt_id",
              default = "C:/Lab/projects/WES/data/bamid_to_bmid_table.txt",
              help = "File containing alternate ID to map sample to phenotype."),
  make_option("--cohort2_alt_id",
              default = NULL,
              help = "File containing alternate ID to map sample to phenotype."),
  make_option("--output_prefix",
              default = file.path(
                "C:/Lab/projects/WES/results/wes_wgs/candidate_gene_filtering",
                "wes_wgs.cadd_phred.maf_0.01"),
              help = "Output directory for GMT parsed files"),
  make_option("--gene_annot",
              default = "C:/Lab/data/genome_wide_gene_list/hg19_gene_coordinates.txt",
              help = "Cutoff for cadd phred score."),
  make_option("--n_resample",
              default = 2000,
              help = "Number of resamplings to perform for empirical p-value"),
  make_option("--genic_intolerance",
              default = "C:/Lab/data/genic_intolerance/RVIS_PLoSGen_EVS_May2013.txt",
              help = "File containing genic intolerance data"),
  make_option("--parallel",
              default = TRUE,
              help = "Logical. Run in parallel."),
  make_option("--cores",
              default = 3,
              help = "If parallel == TRUE, use this number of cores."),
  make_option("--min_count",
              default = c(1, 2),
              help = "Minimum number of affected individuals.")
)

args <- parse_args(OptionParser(option_list = option_list))

# Setup output ----
dir.create(path         = dirname(args$output_prefix),
           showWarnings = FALSE, 
           recursive    = TRUE)

output_prefix <- paste(args$output_prefix, args$gene_set, sep = ".")

cohort1_output_prefix <- dirname(args$output_prefix) %>%
  file.path(basename(tools::file_path_sans_ext(args$cohort1))) %>%
  paste0(".", args$gene_set)

cohort2_output_prefix <- dirname(args$output_prefix) %>%
  file.path(basename(tools::file_path_sans_ext(args$cohort2))) %>%
  paste0(".", args$gene_set)

# Import data ----
message("Importing data")
gene_sets <- read_gmt_file(args$gmt)

message("\tImporting variants")
c1_variants <- read.delim(args$cohort1, header = TRUE, sep = "\t")
c2_variants <- read.delim(args$cohort2, header = TRUE, sep = "\t")

message("\tImporting phenotype")
c1_pheno <- read.table(args$cohort1_pheno, header = TRUE, sep = "\t") %>%
  dplyr::mutate(IID = as.character(IID))
c2_pheno <- read.table(args$cohort2_pheno, header = TRUE, sep = "\t") %>%
  dplyr::mutate(IID = as.character(IID))

message("\tImporting genic intolerance")
genic_intolerance <- import_genic_intolerance(args$genic_intolerance)

message("\tImporting gene annotations")
gene_annotations <- import_gene_annotations(args$gene_annot)

target_gene_set <- gene_sets[[args$gene_set]]
target_gene_set <- target_gene_set[target_gene_set %in% gene_annotations$Gene]

# Extract genes in target gene set ----
message("Identify genes that overlap with target gene set")

# All genes found in the variant file
c1_genes_all <- unique(c1_variants$Gene.refGene)
c2_genes_all <- unique(c2_variants$Gene.refGene)

# Only genes in target gene set
c1_target_genes <- c1_genes_all[c1_genes_all %in% target_gene_set]
c2_target_genes <- c2_genes_all[c2_genes_all %in% target_gene_set]

# Combine the two gene lists
c1_c2_genes <- unique(c(c1_target_genes, c2_target_genes))

# Parse variant annotations ----
message("Extract variants in target gene set")
c1_variants_target_genes <- c1_variants %>%
  dplyr::filter(Gene.refGene %in% c1_target_genes)

c2_variants_target_genes <- c2_variants %>%
  dplyr::filter(Gene.refGene %in% c2_target_genes)

export_table(data     = c1_variants_target_genes,
             filename = paste0(cohort1_output_prefix, ".withcounts"))

export_table(data     = c2_variants_target_genes,
             filename = paste0(cohort2_output_prefix, ".withcounts"))

# Count variant function categories ----
message("Breakdown variants by functional categories")
c1_gene_variant_functions <- count_variant_functions(
  data         = c1_variants, 
  target_genes = c1_target_genes
)

c2_gene_variant_functions <- count_variant_functions(
  data         = c2_variants,
  target_genes = c2_target_genes
)

c1_c2_gene_variant_functions <- c1_gene_variant_functions %>%
  dplyr::full_join(c2_gene_variant_functions, by = "Gene") %>%
  dplyr::rename(Splicing_C1      = Splicing.x,
                Stop_altering_C1 = Stop_altering.x,
                Frameshift_C1    = Frameshift.x,
                Nonsynonymous_C1 = Nonsynonymous.x,
                Total_C1         = Total.x,
                Splicing_C2      = Splicing.y,
                Stop_altering_C2 = Stop_altering.y,
                Frameshift_C2    = Frameshift.y,
                Nonsynonymous_C2 = Nonsynonymous.y,
                Total_C2         = Total.y) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    Splicing_all      = sum(Splicing_C1, Splicing_C2, na.rm = TRUE),
    Stop_altering_all = sum(Stop_altering_C1, Stop_altering_C2, na.rm = TRUE),
    Frameshift_all    = sum(Frameshift_C1, Frameshift_C2, na.rm = TRUE),
    Nonsynonymous_all = sum(Nonsynonymous_C1, Nonsynonymous_C2, na.rm = TRUE),
    Total_all         = sum(Total_C1, Total_C2, na.rm = TRUE)
  )

# Export variant functions
export_table(c1_c2_gene_variant_functions,
             paste(output_prefix, "variant_functions.txt", sep = "."))

# Parse case control for gene set ----
message("Count case-control status for each target gene")
c1_case_control <- count_case_controls(
  data         = c1_variants,
  target_genes = c1_target_genes
)

c2_case_control <- count_case_controls(
  data         = c2_variants,
  target_genes = c2_target_genes
)

# plot all gene case-control counts ----
message("Plot all genes in gene-set")
c1_c2_case_control <- ldply(list(C1 = c1_case_control,
                                 C2 = c2_case_control)) %>%
  dplyr::rename(Cohort = .id)

case_control_melt <- melt(data = dplyr::select(c1_c2_case_control, 
                                               -Total, 
                                               -Variants), 
                          id.vars    = c("Gene", "Cohort"),
                          value.name = "Count") %>%
  dplyr::mutate(Count  = as.integer(Count)) %>%
  dplyr::mutate(Cohort = factor(Cohort, levels = c("C1", "C2"))) %>%
  dplyr::mutate(variable = factor(variable, levels = c("Case", "Control")))

case_control_all_plot <- plot_case_control(
  data     = case_control_melt,
  filename = paste(output_prefix, "all_gene.variant.count.pdf", sep = "."),
  nrow     = NULL,
  width    = 8.5,
  height   = 6,
  units    = "in"
)

# Identify case/control exclusvie gene set ----
message("Extract genes with case/control exclusive counts")
c1_c2_common_genes <- c1_case_control %>%
  dplyr::full_join(c2_case_control, by = "Gene") %>%
  dplyr::left_join(genic_intolerance, by = c("Gene" = "Symbol")) %>%
  dplyr::arrange(RVIS_0.1_percentile) %>%
  dplyr::rename(Case_C1    = Case.x,
                Control_C1 = Control.x,
                Total_C1   = Total.x,
                Case_C2    = Case.y,
                Control_C2 = Control.y,
                Total_C2   = Total.y) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Case    = sum(as.integer(Case_C1), 
                              as.integer(Case_C2),
                              na.rm = TRUE),
                Control = sum(as.integer(Control_C1),
                              as.integer(Control_C2),
                              na.rm = TRUE),
                Total   = sum(as.integer(Total_C1),
                              as.integer(Total_C2),
                              na.rm = TRUE)) %>%
  dplyr::select(Gene, Case_C1, Control_C1, Total_C1, Case_C2, Control_C2,
                 Total_C2, Case, Control, Total, RVIS_0.1, RVIS_0.1_percentile,
                 EdgeCase, OEratio_percentile)

c1_c2_common_genes_two_or_more <- c1_c2_common_genes %>%
  count_exclusive(gene_set = NULL, min_count = 2)

c1_c2_common_genes_one_or_more <- c1_c2_common_genes %>%
  count_exclusive(gene_set = NULL, min_count = 1)

export_table(
  data     = c1_c2_common_genes,
  filename = paste(output_prefix, "combined_gene_count.txt", sep = ".")
)

# Top candidates ----
message("Plot top prioritized genes")
case_control_melt_top <- case_control_melt %>%
  dplyr::filter(Gene %in% c1_c2_common_genes_two_or_more$Gene)
case_control_melt_top[case_control_melt_top == 0] <- NA

top_candidates <- plot_case_control(
  data     = case_control_melt_top,
  filename = paste(output_prefix, "top_gene.variant.count.pdf", sep = "."),
  nrow     = 1,
  width    = 114,
  height   = 50,
  units    = "mm"
)

# Extract variants in top candidates ----
message("Extract variants in top genes")
c1_top_variants <- c1_variants %>%
  dplyr::filter(Gene.refGene %in% c1_c2_common_genes_two_or_more$Gene)

c2_top_variants <- c2_variants %>%
  dplyr::filter(Gene.refGene %in% c1_c2_common_genes_two_or_more$Gene)

c1_one_or_more_variants <- c1_variants %>%
  dplyr::filter(Gene.refGene %in% c1_c2_common_genes_one_or_more$Gene)

c2_one_or_more_variants <- c2_variants %>%
  dplyr::filter(Gene.refGene %in% c1_c2_common_genes_one_or_more$Gene)

export_table(
  data     = c1_top_variants,
  filename = paste(cohort1_output_prefix, "top_variants.withcounts", sep = ".")
)

export_table(
  data     = c2_top_variants,
  filename = paste(cohort2_output_prefix, "top_variants.withcounts", sep = ".")
)

export_table(
  data     = c1_one_or_more_variants,
  filename = paste(cohort1_output_prefix, "one_or_more.withcounts", sep = ".")
)
export_table(
  data     = c2_one_or_more_variants,
  filename = paste(cohort2_output_prefix, "one_or_more.withcounts", sep = ".")
)

# Map to alternate id ----
c1_alt_id <- read.delim(file       = args$cohort1_alt_id, 
                        sep        = "", 
                        header     = TRUE,
                        colClasses = "character")

# Parse individuals that carry the top variants ----
c1_top_individuals <- extract_individuals(
  data      = c1_top_variants,
  phenotype = c1_pheno,
  alt_id    = c1_alt_id
)

c2_top_individuals <- extract_individuals(
  data      = c2_top_variants,
  phenotype = c2_pheno,
  alt_id    = NULL
)

c1_one_or_more_individuals <- extract_individuals(
  data      = c1_one_or_more_variants,
  phenotype = c1_pheno,
  alt_id    = c1_alt_id
)

c2_one_or_more_individuals <- extract_individuals(
  data      = c2_one_or_more_variants,
  phenotype = c2_pheno,
  alt_id    = NULL
)

export_table(
  data     = c1_top_individuals,
  filename = paste(cohort1_output_prefix, "top_individuals.txt", sep = ".")
)
export_table(
  data     = c2_top_individuals,
  filename = paste(cohort2_output_prefix, "top_individuals.txt", sep = ".")
)
export_table(
  data     = c1_one_or_more_individuals,
  filename = paste(cohort1_output_prefix, "one_or_more_individuals.txt", sep = ".")
)
export_table(
  data     = c2_one_or_more_individuals,
  filename = paste(cohort2_output_prefix, "one_or_more_individuals.txt", sep = ".")
)

# Count individuals with variant in gene ----
c1_one_or_more_gene_counts <- count_individuals_per_gene(c1_one_or_more_variants)
c2_one_or_more_gene_counts <- count_individuals_per_gene(c2_one_or_more_variants)

c1_c2_one_or_more_gene_counts <- c1_one_or_more_gene_counts %>%
  dplyr::full_join(c2_one_or_more_gene_counts,
                    by = c("Phenotype" = "Phenotype", "Gene" = "Gene")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Freq = sum(Freq.x, Freq.y, na.rm = TRUE)) %>%
  dplyr::rename(C1_freq = Freq.x, C2_freq = Freq.y)

export_table(data = c1_c2_one_or_more_gene_counts,
             filename = paste(output_prefix,
                              "gene_counts_one_or_more.txt",
                              sep = "."))

# variants and individuals per gene ----
c1_top_variants_select <- c1_top_variants %>%
  dplyr::select(Chr, Start, End, Ref, Alt, Gene.refGene, esp6500siv2_all,
                X1000g2014oct_all, avsnp142, case_id, control_id)
c1_top_variants_select <- adply(c1_top_variants_select, 1, add_novel_variant_id)

c2_top_variants_select <- c2_top_variants %>%
  dplyr::select(Chr, Start, End, Ref, Alt, Gene.refGene, esp6500siv2_all,
                X1000g2014oct_all, avsnp142, case_id, control_id)
c2_top_variants_select <- adply(c2_top_variants_select, 1, add_novel_variant_id)

c1_c2_top_variants_select <- rbind(c1_top_variants_select, 
                                   c2_top_variants_select) %>%
  dplyr::group_by(Gene.refGene) %>%
  dplyr::summarise(variant_id = paste(unique(alt_id), collapse = ";"),
                   case_id = paste(unique(case_id), collapse = ";"),
                   control_id = paste(unique(control_id), collapse = ";")) %>%
  data.frame()

unique_individuals <- extract_all_sample_id(c1_c2_top_variants_select)

# summary tables -----
summary_table <- data.frame(
  Cohort = c("C1", "C2", "Both"),
  Variants = c(nrow(c1_variants_target_genes), 
               nrow(c2_variants_target_genes),
               NA),
  Genes = c(length(c1_target_genes), 
            length(c2_target_genes), 
            length(c1_c2_genes))
)

export_table(
  data     = summary_table,
  filename = paste(output_prefix, "summary_table.txt", sep = ".")
)

################################################################################
# Resampling using all genes to calculate empirical p-value
################################################################################
# Parse case control for all genes ----
c1_variants_parsed <- parse_variant_data(c1_variants)
c2_variants_parsed <- parse_variant_data(c2_variants)

message("Parsing cohort1 dataset")
c1_case_control_all <- count_case_controls(
  data         = c1_variants_parsed,
  target_genes = NULL
)

message("Parsing corhort2 dataset")
c2_case_control_all <- count_case_controls(
  data         = c2_variants_parsed,
  target_genes = NULL
)

# Combine gene counts ----
c1_c2_common_genes_all <- c1_case_control_all %>%
  dplyr::full_join(c2_case_control_all, by = "Gene") %>%
  dplyr::rename(Variants_C1 = Variants.x,
                Case_C1     = Case.x,
                Control_C1  = Control.x,
                Total_C1    = Total.x,
                Variants_C2 = Variants.y,
                Case_C2     = Case.y,
                Control_C2  = Control.y,
                Total_C2    = Total.y) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Case = sum(Case_C1, Case_C2, na.rm = TRUE),
                Control = sum(Control_C1, Control_C2, na.rm = TRUE),
                Total = sum(Total_C1, Total_C2, na.rm = TRUE))

# Run resampling ----
for (min_count in args$min_count){
  resampling_results <- empirical_exclusive(
    target_genes   = target_gene_set,
    all_gene_count = c1_c2_common_genes_all,
    gene_annot     = gene_annotations,
    size_buffer    = 500,
    n_resample     = args$n_resample,
    parallel       = args$parallel,
    cores          = args$cores,
    min_count      = min_count
  )
  
  # pvalue from zscore ----
  resampling_histogram <- plot_resampling(
    data = resampling_results$resampling,
    filename = paste(output_prefix, "min_count", min_count,
                     "histogram.pdf", sep = "."))
  
  zscore_pvalue <- resampling_pvalue(orig        = resampling_results$orig_genes,
                                     resampled   = resampling_results$resampling,
                                     alternative = "one.sided")
  
  resampling_output_files <- export_empirical_pvalue(
    resampling_results = resampling_results,
    zscore_results     = zscore_pvalue,
    output_prefix      = paste(output_prefix, "min_count", min_count, sep = ".")
  )
}
