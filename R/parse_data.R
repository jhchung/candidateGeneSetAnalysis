#' Count the number of cases and controls found for a certain gene.
#' 
#' @param gene \code{character} containing gene name.
#' @param variant_data \code{data.frame} with variant information to parse.
#' @return \code{data.frame} with case/control counts.
#' @export
count_gene_case_control <- function(gene, variant_data){
  temp_vars <- variant_data %>%
    dplyr::filter(Gene.refGene == gene)
  temp_case <- temp_vars %>%
    parse_individuals(id_col = "case_id") %>%
    nrow()
  temp_control <- temp_vars %>%
    parse_individuals(id_col = "control_id") %>%
    nrow()
  temp_total <- temp_case + temp_control
  temp_case_control <- data.frame(Gene     = gene,
                                  Variants = nrow(temp_vars),
                                  Case     = temp_case,
                                  Control  = temp_control,
                                  Total    = temp_total)
  return(temp_case_control)
}

#' Count number of cases and controls.
#' 
#' If \code{target_genes} is \code{NULL}, count all genes. Otherwise, count only
#' genes in \code{target_genes}.
#' 
#' @param data \code{data.frame} containing input data.
#' @param target_genes \code{vector} containing list of genes. default = 
#'   \code{NULL}.
#' @return \code{data.frame} containing case and control count for each gene.
#' @export
count_case_controls <- function(data, target_genes = NULL){
  message("Parsing case-control status")
  
  if(is.null(target_genes)){target_genes <- unique(data$Gene.refGene)}
  
  # Setup count table
  gene_case_control_count <- matrix(nrow = length(target_genes),
                                    ncol = 5) %>%
    data.frame() %>%
    dplyr::rename(Gene     = X1,
                  Variants = X2,
                  Case     = X3,
                  Control  = X4,
                  Total    = X5)
  
  gene_case_control_count <- ldply(.data        = target_genes,
                                   .fun         = count_gene_case_control,
                                   variant_data = data,
                                   .progress    = "text")
  
  gene_case_control_count <- gene_case_control_count %>%
    dplyr::arrange(desc(Total))
  return(gene_case_control_count)
}



#' Count number of individuals who carry variants in a given gene.
#' 
#' @param individuals_to_gene \code{data.frame}.
#' @return \code{data.frame} with individual to gene counts.
#' @export
count_individuals_per_gene <- function(individual_to_gene){
  case_count <- parse_individuals(individual_to_gene, "case_id")
  if (nrow(case_count) > 0) {
    case_count <- case_count %>%
      add_gene_name(individual_to_gene) %>%
      dplyr::select(Gene) %>%
      unlist() %>%
      paste(collapse = ";") %>%
      strsplit(split = ";") %>%
      unlist() %>%
      table() %>%
      data.frame() %>%
      dplyr::rename_(Gene = ".") %>%
      dplyr::mutate(Gene = as.character(Gene))
  } else {
    case_count <- data.frame(Gene = character(), Freq = integer())
  }
  control_count <- parse_individuals(individual_to_gene, "control_id") 
  if (nrow(control_count) > 0){
    control_count <- control_count %>%
      add_gene_name(individual_to_gene) %>%
      dplyr::select(Gene) %>%
      unlist() %>%
      paste(collapse = ";") %>%
      strsplit(split = ";") %>%
      unlist() %>%
      table() %>%
      data.frame() %>%
      dplyr::rename_(Gene = ".") %>%
      dplyr::mutate(Gene = as.character(Gene)) %>%
      dplyr::mutate(Freq = -Freq)
  } else {
    control_count <- data.frame(Gene = character(), Freq = integer())
  }
  all_count <- ldply(list(Case = case_count, Control = control_count)) %>%
    dplyr::rename(Phenotype = .id)
  return(all_count)
}

#' Filter genes that meet case/control exclusive criteria.
#' 
#' @param data \code{data.frame} with gene and subject counts.
#' @param gene_set \code{character} vector with name of genes to check.
#' @param min_count \code{integer} Minimum number of cases or controls in each
#' gene.
#' @return \code{data.frame} with filtered genes.
#' @export
count_exclusive <- function(data, gene_set, min_count = 1){
  require(dplyr)
  if(is.null(gene_set)){gene_set <- data$Gene}
  
  exclusive <- data %>%
    dplyr::filter(Gene %in% gene_set) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Total_case = sum(as.integer(Case_C1),
                                   as.integer(Case_C2),
                                   na.rm = TRUE)) %>%
    dplyr::mutate(Total_control = sum(as.integer(Control_C1),
                                      as.integer(Control_C2),
                                      na.rm = TRUE)) %>%
    dplyr::filter((Total_case >= min_count & Total_control == 0) |
                    (Total_case == 0 & Total_control >= min_count))
  return(exclusive)
}

# test_that("count_exclusive filters genes found in only cases or controls", {
#   expect_equal
# })

#' Extract relevant columns from variant data to reduce data size.
#' 
#' @param variant_data \code{data.frame} containing variant data to process.
#' @return \code{data.frame} containing only relevant columns.
#' @export
parse_variant_data <- function(variant_data){
  parsed_data <- variant_data[, c("Chr", "Start", "End", "Ref", "Alt",
                                  "Func.refGene", "Gene.refGene", 
                                  "GeneDetail.refGene", "ExonicFunc.refGene",
                                  "AAChange.refGene", "case_count", 
                                  "control_count", "case_id", "control_id")]
  return(parsed_data)
}

#' Parse individual IDs
#' 
#' Count the number of times a particular individual shows up in \code{data}.
#' 
#' @param data \code{data.frame} containing input data to parse.
#' @param id_col \code{character} string with name of column to parse.
#' @return \code{data.frame} containing subject ID and the number of times they
#'   are found in \code{data}.
#' @export
parse_individuals <- function(data, id_col){
  data_split <- plyr::alply(data[, id_col], 1, 
                            function(x) strsplit(x, split = ";")[[1]])
  data_split <- as.vector(unlist(data_split))
  
  # Remove "NA"
  data_split <- na.omit(data_split)
  data_split <- data_split[data_split != "NA"]
  
  if (length(data_split) == 0){
    return(data.frame(ID = character(), Count = character()))
  } else {
    data_table <- data.frame(table(data_split))
    names(data_table) <- c("ID", "Count")
    data_table$ID <- as.character(data_table$ID)
    return(data_table)
  }
}

#' Extract all sample IDs from a data.frame.
#' 
#' @inheritParams parse_individuals
#' @return \code{data.frame} with sample ID and number of variants for that
#'   individual carries in \code{data}.
#' @export
extract_all_sample_id <- function(data){
  case_id <- parse_individuals(data, "case_id")
  control_id <- parse_individuals(data, "control_id")
  all_id <- rbind(case_id, control_id)
  return(all_id)
}

#' Extract individuals found in prioritized genes.
#' 
#' @param data \code{data.frame} with variant data.
#' @param phenotype \code{data.frame} containing sample phenotype data.
#' @param alt_id \code{data.frame} with sample alternate ID.
#' @return \code{data.frame} with individual and gene annotations.
#' @export
extract_individuals <- function(data, phenotype, alt_id = NULL){
  individuals <- extract_all_sample_id(data)
  individuals <- dplyr::inner_join(x  = individuals,
                                   y  = phenotype,
                                   by = c("ID" = "IID"))
  if(!is.null(alt_id)){
    individuals <- dplyr::inner_join(x  = individuals,
                                     y  = alt_id,
                                     by = c("ID" = "bamid"))
  }
  individuals <- add_gene_name(individuals = individuals,
                               variants    = data)
  return(individuals)
}

#' Create variant ID if there is no dbSNP rsID
#' 
#' @param data \code{data.frame} with variant information.  
#' @return \code{data.frame} with additional variant ID in \code{alt_id} column.
#' @export
add_novel_variant_id <- function(data){
  if(data$avsnp142 == ""){
    novel_variant <- paste0("chr", data$Chr, ":", data$Start, ":", data$Ref, 
                            ">", data$Alt)
    data$alt_id <- novel_variant
  } else {
    data$alt_id <- data$avsnp142
  }
  return(data)
}

#' Add gene annotations to individuals.
#' 
#' @param individuals \code{data.frame} with individuals to annotate.
#' @param variants \code{data.frame} with variant data.
#' @return \code{data.frame} with additional \code{Gene} column added.
#' @export
add_gene_name <- function(individuals, variants){
  if (nrow(individuals) == 0){
    individuals <- cbind(individuals, Gene = vector(mode = "character"))
    return(individuals)
  } else {
    for (i in 1:nrow(individuals)){
      target_id <- individuals$ID[i]
      gene <- variants[grepl(target_id, variants$case_id), "Gene.refGene"]
      gene <- c(gene, variants[grepl(target_id, variants$control_id), "Gene.refGene"])
      gene <- unique(gene)
      gene <- paste(gene, collapse = ";")
      individuals$Gene[i] <- gene
    }
  }

  return(individuals)
}

get_gene_size_distribution <- function(target_genes, gene_annotations){
  # Process all genes
  all_gene_gene_size <- dplyr::select(gene_annotations, Gene, gene_size)
  all_gene_gene_size <- dplyr::mutate(all_gene_gene_size, Dataset = "Genome")
  
  # Process target gene set
  target_gene_gene_size <- dplyr::filter(gene_annotations, Gene %in% target_genes)
  target_gene_gene_size <- dplyr::select(target_gene_gene_size, Gene, gene_size)
  target_gene_gene_size <- dplyr::mutate(target_gene_gene_size, Dataset = "Target")
  
  # Combine data
  combined_gene_size <- rbind(all_gene_gene_size, target_gene_gene_size) 
  
  gene_size_pvalue <- wilcox.test(x = target_gene_gene_size$gene_size, 
                                  y = all_gene_gene_size$gene_size)
  pvalue_label <- paste("Mann-Whitney p-value =\n", 
                        signif(gene_size_pvalue$p.value, digits = 3))
  
  gene_size_results <- list(pvalue        = gene_size_pvalue,
                            pvalue_label  = pvalue_label,
                            gene_size     = combined_gene_size)
  class(gene_size_results) <- "gene_size_data"
  return(gene_size_results)
}