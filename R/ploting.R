#' Create violin plots of gene size across the whole genome vs in the target
#' gene set.
#' 
#' Takes the output from \code{\link{get_gene_size_distribution}}.
#' 
#' @param gene_size_data \code{gene_size_data} object from 
#'   \code{\link{get_gene_size_distribution}}
#' @param output_prefix \code{character} output file prefix.
#' @return ggplot2 plot
#' @seealso \code{\link{get_gene_size_distribution}}.
#' @export
plot_gene_size <- function(gene_size_data, output_prefix){
  if(class(gene_size_data) != "gene_size_data"){
    stop("Wrong input data class: gene_size_data\nPlease use output from get_gene_size_distribution")
  }
  
  output_file <- paste(output_prefix, "gene_size.pdf", sep = ".")
  
  gene_size_plot <- ggplot(data = gene_size_data$gene_size, 
                           aes(factor(Dataset), log10(gene_size))) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.1) +
    ggplot2::annotate("text", x = 1.5, y = 2, 
                      label = gene_size_data$pvalue_label) + 
    ggplot2::ggtitle("Comparison of gene sizes") +
    ggplot2::xlab("Gene set") +
    ggplot2::ylab("log10(Gene size [bp])") +
    ggplot2::ggsave(filename = output_file,
                    width = 183, height = 100, units = "mm")
  return(gene_size_plot)
}

#' Plot case control bar graphs.
#' 
#' @param data \code{data.frame} formatted case/control data.
#' @inheritParams ggplot2::facet_wrap
#' @inheritParams ggplot2::ggsave
#' @export
plot_case_control <- function(data, filename, nrow = NULL, width = 8.5, 
                              height = 6, units = "in"){
  purple_colors <- c("#BDBAD8", "#6258B0")
  grey_colors   <- c("#7A7A7A", "#C9C9C9")
  green_red     <- c("#FF9999", "#99FF99")
  case_control_colors <- green_red
  
  case_control_plot <- 
    ggplot2::ggplot(data, aes(x = Cohort, y = Count, fill = factor(variable))) +
    ggplot2::geom_bar(stat='identity', position = position_dodge()) +
    ggplot2::scale_fill_manual(name = "Phenotype",
                               values = case_control_colors) + 
    ggplot2::facet_wrap(~ Gene, nrow = nrow) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   strip.background = element_blank(),
                   strip.text       = element_text(size = 9),
                   axis.text        = element_text(size = 7),
                   axis.title       = element_text(size = 9),
                   legend.text      = element_text(size = 9),
                   legend.title     = element_text(size = 9),
                   panel.border     = element_rect(color = "grey"),
                   legend.position  = "bottom",
                   legend.key.size  = grid::unit(3, "mm"),
                   legend.margin    = grid::unit(c(-36, 0, 0, 0), "mm"),
                   plot.margin      = grid::unit(c(1, 1, -3, 1), "mm")) +
    ggplot2::ggsave(filename = filename,
           width = width, height = height, units = units)
  return(case_control_plot)
}

#' Plot histogram of resampling results.
#' 
#' @param data \code{vector} of resampling results.
#' @inheritParams ggplot2::ggsave
#' @return ggplot object
#' @export
plot_resampling <- function(data, filename){
  data_df <- data.frame(n_genes = data)
  resampling_histogram <- ggplot2::ggplot(data = data_df, aes(x = n_genes)) +
    ggplot2::geom_histogram(aes(y = ..density..), binwidth = 1) + 
    ggplot2::theme_bw()
  ggsave(filename, width = 183, height = 100, units = "mm")
  return(resampling_histogram)
}