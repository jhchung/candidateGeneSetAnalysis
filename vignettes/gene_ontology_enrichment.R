# Gene Ontology enrichment analysis

# Install bioconductor packages ----
# source("http://bioconductor.org/biocLite.R")
# biocLite("topGO")
# biocLite("org.Hs.eg.db")
# biocLite("Rgraphviz")

# Setup ----
options(stringsAsFactors = FALSE)
library(topGO)
library(dplyr)
library(plyr)
library(org.Hs.eg.db)
library("Rgraphviz")

# Functions -----

#' setup output file.
create_output_file <- function(input_file, output_dir, pattern, replacement){
  output_file <- basename(input_file)
  output_file <- gsub(pattern = pattern, replacement = replacement, output_file)
  output_file <- file.path(output_dir, output_file)
  return(output_file)
}

create_target_gene_list <- function(target_genes, background_genes){
  target_gene_list <- factor(as.integer(background_genes %in% target_genes))
  names(target_gene_list) <- background_genes
  return(target_gene_list)
}

setup_topGOdata_objects <- function(gene_list, ontology = c("BP", "MF", "CC"), 
                                    node_size = 10, description = "GO analysis",
                                    mapping = "org.Hs.eg.db", ID = "Symbol"){
  go_data_object <- new(
    Class = "topGOdata",
    ontology = ontology,
    allGenes = gene_list,
    nodeSize = node_size,
    description = description,
    annotationFun = annFUN.org,
    mapping = mapping,
    ID = ID
  )
  return(go_data_object)
}

extract_tables <- function(go_data_object_list, go_results_list){
  go_tables <- list()
  for (target_category in names(go_data_object_list)){
    message("Extracting tables for: ", target_category)
    go_tables[[target_category]] <- topGO::GenTable(
      object   = go_data_object_list[[target_category]],
      Fisher   = go_results_list[[target_category]],
      topNodes = length(go_results_list[[target_category]]@score),
      numChar  = 1000
    )
  }
  return(go_tables)
}

add_gene_names <- function(result_table_list, go_data_list, gene_lists){
  annotated_tables <- list()
  for (target_genes in names(gene_lists)){
    message("Adding gene names for: ", target_genes)
    annotated_tables[[target_genes]] <- get_significant_genes(
      go_object    = go_data_list[[target_genes]],
      result_table = result_table_list[[target_genes]],
      gene_set     = gene_lists[[target_genes]])
  }
  return(annotated_tables)
}

#' Wrapper to run topGO using multiple term weighting algorithms
#' 
#' @args go_object \code{topGOdata} object. 
#' @return named list of \code{topGOresult} objects.
run_multiple_test <- function(go_object){
  message("Running classic test")
  classic_test <- topGO::runTest(
    object    = go_object,
    algorithm = "classic",
    statistic = "fisher")
  
  message("Running elim test")
  elim_test <- topGO::runTest(
    object    = go_object,
    algorithm = "elim",
    statistic = "fisher")
  
  message("Running weight test")
  weight_test <- topGO::runTest(
    object    = go_object,
    algorithm = "weight",
    statistic = "fisher")
  
  message("Running weight01 test")
  weight01_test <- topGO::runTest(
    object    = go_object,
    algorithm = "weight01",
    statistic = "fisher")
  
  message("Running parentchild test")
  parentchild_test <- topGO::runTest(
    object    = go_object,
    algorithm = "parentchild",
    statistic = "fisher")
  
  return(list(classic     = classic_test,
              elim        = elim_test,
              weight      = weight_test,
              weight01    = weight01_test,
              parentchild = parentchild_test))
}

get_significant_genes <- function(go_object, result_table, gene_set){
  significant_genes <- c()
  for (i in 1:nrow(result_table)){
    go_term_genes <- topGO::genesInTerm(go_object, result_table[i, "GO.ID"])[[1]]
    go_gene_set_intersect <- go_term_genes[go_term_genes %in% gene_set]
    go_gene_set_intersect <- paste(go_gene_set_intersect, collapse = ";")
    significant_genes <- c(significant_genes, go_gene_set_intersect)
  }
  output_table <- result_table %>%
    dplyr::mutate(Significant_genes = significant_genes)
  return(output_table)
}

export_top_go <- function(result_table, output_file){
  require(dplyr)
  output_table <- result_table %>%
    dplyr::rename(Prioritized       = Significant,
                  Prioritized_genes = Significant_genes)
  write.table(x         = output_table,
              file      = output_file,
              sep       = "\t",
              quote     = FALSE,
              row.names = FALSE,
              col.names = TRUE)
}

#' Combine pvalues from different algorithms as the mean of pvalues.
#' 
#' @args pvalues \code{vector} of pvalues to comibine.
#' @return combined pvalues
combine_pvalues <- function(pvalues){
  combined_pvalue <- exp(1 / length(pvalues) * sum(log(as.numeric(pvalues))))
  return(combined_pvalue)
}


#' Adjust pvalues for multiple testing.
#' 
#' @args data \code{data.frame} with data to adjust.
#' @args method \code{character} string containing the method to adjust p-values. See \code{p.adjust} for details.
#' @args cutoff \code{numeric}. Adjusted p-value cutoff.
#' @return \code{data.frame} with additional column containing adjusted p-value.
adjust_multiple_testing <- function(data, method = "BY", cutoff = 0.05){
  adjusted_data <- data
  adjusted_data$adjusted_pvalue <- stats::p.adjust(data$Fisher, method = method)
  if (is.na(cutoff)){
    return(adjusted_data)
  } else {
    adjusted_data <- base::subset(adjusted_data, subset = adjusted_pvalue < cutoff)
    return(adjusted_data)
  }
}

export_go_tables <- function(go_tables, input_file, output_dir, pattern, go){
  for (target_set in names(go_tables)){
    output_pattern <- paste0(target_set, "_go_", go, "_annotated.txt")
    output_file <- create_output_file(
      input_file  = input_file,
      output_dir  = output_dir,
      pattern     = pattern,
      replacement = output_pattern
    )
    message("Export data to: ", output_file)
    export_top_go(result_table = go_tables[[target_set]],
                  output_file = output_file)
  }
}

# Input data ----
gene_count_file <- "results/wes_wgs/candidate_gene_filtering/wes_wgs.cadd_phred.maf_0.01.tbx1_pathway_animal_model.combined_gene_count.txt"
wes_variant_file <- "results/wes186/variant_filtering/22q_exomes.cadd_phred.maf_0.01.withcounts"
wgs_variant_file <- "results/wgs100/variant_filtering/22q_batch1.cadd_phred.maf_0.01.withcounts"
pathway_file <- "data/gene_lists/tbx1_pathway_from_animal_models_human_homolog.txt"

# Parameters ----
output_dir <- "results/wes_wgs/gene_ontology"
input_file_pattern <- "combined_gene_count.txt" 
algorithm <- "weight01"

go_bp <- TRUE
go_cc <- FALSE
go_mf <- TRUE
ontology_list <- c("CC", "MF", "BP")

# Setup output ----
output_dir <- file.path(output_dir,
                        paste0("algorithm_", algorithm))
dir.create(output_dir)

# Import ----
gene_count    <- read.table(gene_count_file,  header = TRUE, sep = "\t")
wes_variant   <- read.delim(wes_variant_file, header = TRUE, sep = "\t")
wgs_variant   <- read.delim(wgs_variant_file, header = TRUE, sep = "\t")
pathway_genes <- read.table(pathway_file,     header = TRUE, sep = "\t")

# Setup gene universe ----
wes_genes <- wes_variant %>%
  dplyr::select(Gene.refGene) %>%
  unlist() %>%
  unique()
wgs_genes <- wgs_variant %>%
  dplyr::select(Gene.refGene) %>%
  unlist() %>%
  unique()
all_genes <- unique(c(wes_genes, wgs_genes))

# Number of unique variants ----
wes_variant_positions <- wes_variant %>%
  dplyr::select(Chr, Start, End, Ref, Alt) %>%
  plyr::alply(.margins = 1, .fun = paste, collapse = ":") %>%
  unlist()
wgs_variant_positions <- wgs_variant %>%
  dplyr::select(Chr, Start, End, Ref, Alt) %>%
  plyr::alply(.margins = 1, .fun = paste, collapse = ":") %>%
  unlist()

wes_wgs_variant_positions <- c(wes_variant_positions, wgs_variant_positions) %>%
  unique()

# filter genes ----
control_only <- gene_count %>%
  dplyr::rowwise() %>%
  dplyr::filter(Control >= 1 & Case == 0)
case_only <- gene_count %>%
  dplyr::rowwise() %>%
  dplyr::filter(Control == 0 & Case >= 1)

# target gene list ----
target_genes <- list(
  Case         = case_only$Gene,
  Control      = control_only$Gene,
  Case_Control = c(case_only$Gene, control_only$Gene),
  All          = gene_count$Gene
)

gene_lists <- llply(.data            = target_genes,
                    .fun             = create_target_gene_list,
                    background_genes = all_genes)

# Run topGO test ----
# Biological Process
for (ontology in ontology_list){
  message("Running enrichment analysis for ontology: ", ontology)
  go_data_objects <- llply(gene_lists, setup_topGOdata_objects,
                              ontology = ontology, node_size = 10)
  
  if (algorithm == "multiple"){
    go_results <- llply(.data = go_data_objects,
                           .fun  = run_multiple_test)
    go_tables <- list()
    for (target_list in names(go_results)){
      message("Collecting tables for: ", target_list)
      go_tables[[target_list]] <- topGO::GenTable(
        object      = go_data_objects[[target_list]],
        classic     = go_results[[target_list]]$classic,
        elim        = go_results[[target_list]]$elim,
        weight      = go_results[[target_list]]$weight,
        weight01    = go_results[[target_list]]$weight01,
        parentchild = go_results[[target_list]]$parentchild,
        ranksOf     = "classic",
        topNodes    = length(go_results[[target_list]]$classic@score),
        numChar     = 1000
      )
      
      go_tables[[target_list]] <- go_tables[[target_list]] %>%
        dplyr::rowwise() %>%
        dplyr::mutate(combined_pvalue = combine_pvalues(
          c(classic, elim, weight, parentchild)
        )) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(adjusted_pvalue = p.adjust(p      = combined_pvalue, 
                                                 method = "BY")) %>%
        dplyr::arrange(adjusted_pvalue)
      
      go_tables[[target_list]] <- get_significant_genes(
        go_object    = go_data_objects[[target_list]],
        result_table = go_tables[[target_list]],
        gene_set     = target_genes[[target_list]]
      )
    }
    export_go_tables(go_tables  = go_tables,
                     input_file = gene_count_file,
                     output_dir = output_dir,
                     pattern    = input_file_pattern,
                     go         = ontology)
  } else {
    go_results <- llply(.data      = go_data_objects,
                        .fun       = topGO::runTest,
                        algorithm  = algorithm,
                        statistic  = "fisher")
    go_tables <- extract_tables(go_data_object_list = go_data_objects,
                                go_results_list     = go_results)
    go_tables_annotated <- llply(.data  = go_tables,
                                 .fun   = adjust_multiple_testing,
                                 method = "BY",
                                 cutoff = NA)
    go_tables_annotated <- add_gene_names(
      result_table_list = go_tables_annotated,
      go_data_list      = go_data_objects,
      gene_lists        = target_genes
    )
    export_go_tables(go_tables  = go_tables_annotated,
                     input_file = gene_count_file,
                     output_dir = output_dir,
                     pattern    = input_file_pattern,
                     go         = ontology)
  }
}

# # Plot GO graph ----
# topGO::printGraph(
#   object        = case_go_bp,
#   result        = case_bp_results,
#   firstSigNodes = nrow(case_bp_results_table_fdr), 
#   useInfo       = "all",
#   fn.prefix     = file.path(output_dir, "GOgraph_case_bp"), 
#   pdfSW         = TRUE
# )
# 
# topGO::printGraph(
#   object        = control_go_bp,
#   result        = control_bp_results,
#   firstSigNodes = 5, 
#   useInfo       = "all",
#   fn.prefix     = file.path(output_dir, "GOgraph_control_bp"), 
#   pdfSW         = TRUE
# )
# 
# topGO::printGraph(
#   object        = case_control_go_bp,
#   result        = case_control_bp_results,
#   firstSigNodes = 5, 
#   useInfo       = "all",
#   fn.prefix     = file.path(output_dir, "GOgraph_case_control_bp"), 
#   pdfSW         = TRUE
# )
