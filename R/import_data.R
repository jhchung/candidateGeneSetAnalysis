#' Import and format gene annotations.
#' 
#' @param annotation_file \code{character}. Path to gene annotation file.
#' @return \code{data.frame} with gene annotations.
#' @export
import_gene_annotations <- function(annotation_file){
  gene_annotations <- read.delim(annotation_file, sep = "\t", header = TRUE) %>%
    dplyr::rename(Gene  = hg19.kgXref.geneSymbol,
                  chrom = hg19.knownGene.chrom,
                  start = hg19.knownGene.txStart,
                  end   = hg19.knownGene.txEnd) %>%
    dplyr::mutate(center = (start + end / 2) / 1e6)  %>%
    dplyr::group_by(Gene) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(chrom = as.integer(gsub("chr", "", chrom))) %>%
    dplyr::mutate(gene_size = end - start)
  return(gene_annotations)
}

#' Import genic intolerance data into a \code{data.frame}.
#' 
#' Data downloaded from http://igm.cumc.columbia.edu/GenicIntolerance/data/SCORES_n12_4NR_v16May15.txt
#' 
#' @param input_file \code{character}. Path to genic intolerance file
#' @return \code{data.frame} containing RVIS data.
#' @export
import_genic_intolerance <- function(input_file){
  data <- read.table(file   = input_file,
                     sep    = "\t",
                     header = TRUE) %>%
    dplyr::rename(Symbol              = CCDS_r9,
                  RVIS_0.1            = X0.1.RVIS,
                  RVIS_0.1_percentile = X0.1.RVIS.,
                  OEratio_percentile  = OEratio.tile) %>%
    dplyr::select(Symbol, 
                  RVIS_0.1,
                  RVIS_0.1_percentile,
                  EdgeCase, 
                  OEratio_percentile)
  return(data)
}

#' Parse input \code{data.frame} for variant function and count the total in each category.
#' 
#' @param data \code{data.frame} containing variant information.
#' @return \code{data.frame} with variant function and count.
#' @export
parse_variant_functions <- function(data){
  splicing      <- 0
  stop_altering <- 0
  frameshift    <- 0
  nonsynonymous <- 0
  
  for (i in 1:nrow(data)){
    if (grepl("splicing", data$Func.refGene[[i]])){
      splicing <- splicing + 1
    } else if (grepl("stop", data$ExonicFunc.refGene[[i]])){
      stop_altering <- stop_altering + 1
    } else if (grepl("frameshift", data$ExonicFunc.refGene[[i]])){
      frameshift <- frameshift + 1
    } else if (grepl("nonsynonymous", data$ExonicFunc.refGene[[i]])){
      nonsynonymous <- nonsynonymous + 1
    }
  }
  total <- sum(splicing, 
               stop_altering, 
               frameshift, 
               nonsynonymous)
  return(data.frame(Splicing      = splicing,
                    Stop_altering = stop_altering,
                    Frameshift    = frameshift,
                    Nonsynonymous = nonsynonymous,
                    Total         = total))
}

#' For a given set of genes, summarize the number of variants and their 
#' functions.
#' 
#' @param data \code{data.frame} containing input data from 
#'   \code{\link{parse_variant_functions}}.
#' @param target_genes \code{vector} containing target genes to summarize.
#' @return \code{data.frame} continaing summary of variant functions for each 
#'   gene.
#' @seealso \code{\link{parse_variant_functions}}
#' @export
count_variant_functions <- function(data, target_genes){
  gene_variant_functions <- data.frame()
  for (i in 1:length(target_genes)){
    temp_vars <- dplyr::filter(data, Gene.refGene == target_genes[[i]])
    temp_categories <- parse_variant_functions(temp_vars)
    temp_categories <- cbind(Gene = target_genes[[i]],
                             temp_categories)
    gene_variant_functions <- rbind(gene_variant_functions,
                                    temp_categories)
  }
  
  return(dplyr::arrange(gene_variant_functions, desc(Total)))    
}


#' Format \code{.gmt} data for use in \code{SNPath}
#' 
#' @param gmt_data \code{.gmt} data as read using \code{scan}.
format_gmt <- function(gmt_data){
  require(plyr)
  formatted_gmt_data <- llply(gmt_data, function(x) as.array(x[3:length(x)]))
  set_names <- laply(gmt_data, function(x) x[[1]])
  
  names(formatted_gmt_data) <- set_names
  return(formatted_gmt_data)
}


#' Import GMT formatted file for use in \code{SNPath}.
#' 
#' @param gmt_file Path to input \code{.gmt} file.
#' @export
read_gmt_file <- function(gmt_file){
  message("Reading .gmt file")
  gmt_data <- scan(gmt_file, what = "character", sep = "\n", fill = TRUE, quote = "")
  gmt_data <- strsplit(gmt_data, "\t")
  gmt_data <- format_gmt(gmt_data)
  return(gmt_data)
}

