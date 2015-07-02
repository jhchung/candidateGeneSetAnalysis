#' Filter previously selected genes from input \code{data.frame}.
#' 
#' @param input \code{data.frame} with genes to check.
#' @param previous \code{data.frame} with previously selected genes.
#' @return \code{data.frame} with filtered results.
remove_prevously_selected_genes <- function(input, previous){
  if (nrow(previous) == 0){
    return(input)
  } else {
    selected_genes <- previous$Gene
    filtered_data <- dplyr::filter(input, !Gene %in% selected_genes)
    return(filtered_data)
  }
}

#' Random set of genes with similar size as those in \code{target_genes}.
#' 
#' Select random subset of rows where nrows = number of target genes without
#' replacement.
#' 
#' @param target_genes \code{vector} containing gene names in target set.
#' @param background_genes \code{data.frame} containing all genes to search.
#' @param size_buffer \code{integer}. Buffer to pick genes of target size in bp.
#' @export
get_random_set <- function(target_genes, background_genes, size_buffer = 200){
#   # Filter background genes 
#   background_no_target <- dplyr::filter(background_genes, 
#                                         !Gene %in% target_genes)
  random_subset <- data.frame()
  for (i in 1:length(target_genes)){
    if(is.na(size_buffer)){
      # Select any random subset of rows 
      random_subset <- background_genes[sample(1:nrow(background_genes), 
                                               length(target_genes)), ]
    } else {
      target_gene_annotation <- dplyr::filter(background_genes,
                                              Gene == target_genes[i])
      min_gene_size <- target_gene_annotation$gene_size - size_buffer
      max_gene_size <- target_gene_annotation$gene_size + size_buffer
      
      # Filter based on gene size.
      size_subset <- dplyr::filter(
        background_genes,
        gene_size >= min_gene_size & gene_size <= max_gene_size
      )
      
      # Keep track of genes already selected
      size_subset <- remove_prevously_selected_genes(size_subset, random_subset)

      random_gene <- size_subset[sample(1:nrow(size_subset), 1),]
      random_subset <- rbind(random_subset, random_gene)
    }
  }
  return(random_subset)
}

#' Get count of random gene sets that pass filtering criteria.
#' 
#' @param target_genes \code{vector} containing gene names in target set.
#' @param all_gene_count \code{data.frame} with count of samples for all genes.
#' @param gene_annot \code{data.frame with gene annotations}.
#' @param size_buffer \code{integer} buffer in base pairs to search for genes of
#'   similar size to input gene.
#' @param n_resample \code{integer} number of resamplings to perform.
#' @param parallel \code{logical} run analysis in parallel?
#' @param min_count \code{integer} minimum number of individuals affected for a
#'   given gene.
#' @return \code{vector} with number of genes that pass filtering criteria.
#' @export
resample_counts <- function(target_genes, all_gene_count, gene_annot, 
                            size_buffer, n_resample = 500, 
                            parallel = FALSE, cores = 1, min_count = 1){
  require(dplyr)
  if (parallel){
    # Setup parallel environment
    snow_cl <- snow::makeCluster(cores)
    snow::clusterExport(snow_cl, c("get_random_set", "count_exclusive"))
    doSNOW::registerDoSNOW(snow_cl)
    
    # Start timer
    message("Resampling started!")
    message("Monitor progress at ", file.path(getwd(), "resampling_log.txt"))
    cat(paste(Sys.time(), "Starting resampling!!\n"), file = "resampling_log.txt")
    
    all_resampled <- foreach::foreach(i = 1:n_resample, .combine = "c") %dopar% {                             
      
      if ( i %% 50 == 0){
        cat(paste(Sys.time(), "Running iteration:", i, "of", n_resample, "\n"),
            file = "resampling_log.txt", append = TRUE)
      }
      
      random_subset <- get_random_set(
        target_genes     = target_genes,
        background_genes = gene_annot,
        size_buffer      = size_buffer
      )
      
      random_set_exclusive <- count_exclusive(
        data      = all_gene_count,
        gene_set  = random_subset$Gene,
        min_count = min_count
      )
      
      nrow(random_set_exclusive)
    }
    
    # Stop cluster
    snow::stopCluster(snow_cl)
    message("Resampling finished!!")
    cat(paste(Sys.time(), "Finished!!\n"), file = "resampling_log.txt", 
        append = TRUE)
    
  } else {
    all_resampled <- rep(NA, n_resample)
    pb <- txtProgressBar(min = 0, max = n_resample, style = 3)
    for(i in 1:n_resample){
      setTxtProgressBar(pb, i)
      random_subset <- get_random_set(target_genes     = target_genes,
                                      background_genes = gene_annot,
                                      size_buffer      = size_buffer)
      random_set_exclusive <- count_exclusive(data      = all_gene_count,
                                              gene_set  = random_subset$Gene,
                                              min_count = min_count)
      all_resampled[i] <- nrow(random_set_exclusive)
    }
    close(pb)
  }
  return(all_resampled)
}

#' Calculate empirical p-value from resampled data.
#' 
#' @param orig \code{integer}. Number of genes found exclusively in cases or
#'   controls.
#' @param resampled \code{vector} containing number of genes found exclusively
#'   in cases or controls in a random set of genes.
#' @return \code{numeric}. The empirical p-value
calculate_empirical_pvalue <- function(orig, resampled){
  more_than_orig   <- resampled[resampled >= orig]
  more_than_orig   <- length(more_than_orig)
  empirical_pvalue <- (more_than_orig + 1) / (length(resampled) + 1)
  return(empirical_pvalue)
}

#' Calculate empirical p-value for number of random gene sets with the same or
#' more "prioritized" genes than the target gene set.
#' 
#' @param target_genes \code{data.frame} containing gene names in target gene
#'   set.
#' @param all_gene_count \code{data.frame} of sample counts for all genes.
#' @param gene_annot \code{data.frame} of gene annotations.
#' @param size_buffer \code{integer} gene size buffer in base pairs to select a
#'   random set of genes.
#' @param size_buffer \code{integer} buffer in base pairs to search for genes of
#'   similar size to input gene.
#' @param n_resample \code{integer} number of resamplings to perform.
#' @param parallel \code{logical} run analysis in parallel?
#' @param cores \code{integer} number of cores to use if \code{parallel == TRUE}.
#' @param min_count \code{integer} minimum number of individuals affected for a
#'   given gene.
#' @return \code{list} containing empirical p-value and all resampling results.
#' @export
empirical_exclusive <- function(target_genes, all_gene_count, gene_annot,
                                size_buffer, n_resample, parallel = FALSE,
                                cores = 1, min_count = 1){
  message("Min count: ", min_count)
  target_exclusive <- count_exclusive(data      = all_gene_count,
                                      gene_set  = target_genes,
                                      min_count = min_count)
  n_genes_orig <- nrow(target_exclusive)
  
  all_resampled <- resample_counts(target_genes   = target_genes,
                                   all_gene_count = all_gene_count,
                                   gene_annot     = gene_annot,
                                   size_buffer    = size_buffer,
                                   n_resample     = n_resample,
                                   parallel       = parallel,
                                   cores          = cores,
                                   min_count      = min_count)
  
  empirical_pvalue <- calculate_empirical_pvalue(orig      = n_genes_orig,
                                                 resampled = all_resampled)
  return(list(pvalue     = empirical_pvalue,
              orig_genes = n_genes_orig,
              resampling = all_resampled))
}

#' Calculate p-value using z-scores.
#' 
#' Assuming the resampled results form a normal distribution. Use z-scores to
#' calculate p-value of the original result.
#' 
#' @param orig \code{integer} number of genes in the original gene set with 
#'   variants found exclusively in cases or controls.
#' @param resampled \code{vector} of integers where each element in the vector
#'   is a random gene set the number of genes with
#'   variants found exclusively in controls.
#' @return \code{numeric} pvalue of \code{orig}.
#' @export
resampling_pvalue <- function(orig, resampled, 
                              alternative = c("one.sided", "two.sided")){
  all_zscore <- base::scale(c(orig, resampled), center = TRUE, scale = TRUE)
  orig_zscore <- all_zscore[1]
  if (alternative == "one.sided"){
    pvalue <- pnorm(-abs(orig_zscore))
  } else if (alternative == "two.sided"){
    pvalue <- 2 * pnorm(-abs(orig_zscore))
  }
  return(pvalue)
}

