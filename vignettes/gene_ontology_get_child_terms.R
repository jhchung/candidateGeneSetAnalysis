# Extract genes in GO terms

# Setup ----
library("org.Hs.eg.db")
library("GO.db")
library(annotate)
library(plyr)
library(dplyr)

# Functions ----
#' Get vector of additional GO terms
get_extra_terms <- function(go_term, which_terms = "none"){
  term_ontology <- GO.db::GOTERM[[go_term]]@Ontology
  if(which_terms == "children"){
    if (term_ontology == "BP"){
      extra_terms  <- GO.db::GOBPCHILDREN[[go_term]]
    } else if (term_ontology == "MF") {
      extra_terms  <- GO.db::GOMFCHILDREN[[go_term]]
    } else if (term_ontology == "CC"){
      extra_terms  <- GO.db::GOCCCHILDREN[[go_term]]
    }
  } else if (which_terms == "offspring"){
    if (term_ontology == "BP"){
      extra_terms <- GO.db::GOBPOFFSPRING[[go_term]]
    } else if (term_ontology == "MF") {
      extra_terms <- GO.db::GOMFOFFSPRING[[go_term]]
    } else if (term_ontology == "CC"){
      extra_terms <- GO.db::GOCCOFFSPRING[[go_term]]
    }
  } else {
    extra_terms = NULL
  }
  all_terms <- unique(c(go_term, extra_terms))
  return(all_terms)
}

collect_genes <- function(go_term, go_object, evidence = "all"){
  genes <- go_object[[go_term]]
  if (evidence != "all"){
    genes <- genes[names(genes) == evidence]
  }
  genes <- unique(as.vector(genes))
  genes <- list(genes)
  names(genes) <- go_term
  return(genes)
}

export_gmt <- function(go_term, genes, extra_terms, out_file, append = FALSE){
  term_info <- GO.db::GOTERM[[go_term]]
  if (extra_terms != "none"){
    gene_set_name <- paste(go_term, "and", extra_terms, sep = "_")
    gene_set_name <- gsub(":", "_", gene_set_name)
  } else {
    gene_set_name <- gsub(":", "_", go_term)
  }
  description <- term_info@Term
  gmt_output <- c(gene_set_name, description, genes)
  gmt_output <- paste(gmt_output, collapse = "\t")
  write(gmt_output, sep = "\n", file = out_file, append = append)
}

# Parameters ----
go_term <- "GO:2000649"
go_name <- GO.db::GOTERM[[go_term]]@Term %>%
  gsub(pattern = " ", replacement = "_")

# One of c("offspring", "children", "none")
which_extra_terms <- c("offspring", "children", "none")

output_path <- "results/wes_wgs/gene_ontology_gene_set"

# Create output file ----
output_file <- file.path(output_path, paste0(go_name, ".gmt"))
dir.create(dirname(output_file), recursive = TRUE)
file.create(output_file)

# Search gene ontology ----
for (extra_term in which_extra_terms){
  # Collect child terms ----
  all_go_terms <- get_extra_terms(go_term, extra_term)
  
  # Setup GO objects ----
  # setup entrez object
  entrez_object <- org.Hs.egGO
  
  # Map GO terms to Entrez gene ID ----
  go_object <- AnnotationDbi::as.list(org.Hs.egGO2EG)
  
  
  all_term_entrez <- laply(all_go_terms, collect_genes, go_object = go_object,
                           evidence = "all")
  all_term_entrez <- all_term_entrez %>%
    as.vector() %>%
    unlist() %>%
    unique()
  
  all_term_symbol <- annotate::getSYMBOL(all_term_entrez, data = "org.Hs.eg") %>%
    unlist() %>%
    unique()
  
  
  export_gmt(go_term, genes = all_term_symbol, extra_terms = extra_term, 
             out_file = output_file, append = TRUE)
  
}
