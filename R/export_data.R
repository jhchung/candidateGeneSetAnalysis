#' Helper function to export tables.
#' 
#' @param data \code{data.frame} with data to export.
#' @param filename \code{character} string containing output file path.
#' @export
export_table <- function(data, filename){
  write.table(x = data,
              file = filename,
              quote = FALSE,
              col.names = TRUE,
              row.names = FALSE,
              sep = "\t")
}

#' Export resampling results.
#' 
#' @param resampling_results \code{list} with empirical p-value and resampling
#'   values.
#' @param results from p-value calculations using z-scores.
#' @param output_prefix \code{character} inditcating output file prefix.
#' @return \code{list} with output file names for pvalue and resampling data.
#' @export
export_empirical_pvalue <- function(resampling_results, zscore_results,
                                    output_prefix){
  empirical_pvalue_file <- paste(output_prefix,  "empirical_pvalue.txt", 
                                 sep = ".")
  resampling_file <- paste(output_prefix, "resampling.txt", sep = ".")
  
  write(paste("empirical p-value =", resampling_results$pvalue),
        file = empirical_pvalue_file)
  write(paste("z score p-value =", zscore_results),
        file = empirical_pvalue_file, append = TRUE)
  write(paste("Number of prioritized genes =", resampling_results$orig_genes),
        file = empirical_pvalue_file, append = TRUE)
  write(paste("n-resampling =", length(resampling_results$resampling)),
        file = empirical_pvalue_file, append = TRUE)
  
  write(paste("# Resampling values:\n",
              paste(resampling_results$resampling, collapse = "\t"), sep = ""),
        file = resampling_file)
  return(list(empirical_pvalue = empirical_pvalue_file,
              resampling       = resampling_file))
}