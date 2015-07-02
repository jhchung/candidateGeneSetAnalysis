# Plot average read depth in 22q11.2 deleted region

library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
# source("http://bioconductor.org/workflows.R")
# workflowInstall("liftOver")
library(rtracklayer)
library(GenomicRanges)
library(AnnotationHub)
library(extrafont)
font_import(paths = "C:/Windows/Fonts", pattern = "arial", prompt = FALSE)
loadfonts()

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

segdup_color <- function(data){
  if(data >= 0.9 & data < 0.98){
    segment_group <- 1
  } else if (data >= 0.98 & data < 0.99){
    segment_group <- 2
  } else if (data >= 0.99){
    segment_group <- 3
  }
  return(segment_group)
}

################################################################################
wes_file <- "results/wes186/22q_region/22q_exomes.22q_region.ldepth.mean"
wgs_file <- "results/wgs100/22q_region/22q_batch1.22q_region.ldepth.mean"
segdup_file <- "data/segmental_duplications/chr22_segmental_dups.hg19.txt"

output_file <- "results/wes_wgs/22q_region/wes_wgs.22q_region.ldepth.mean.pdf"

dir.create(dirname(output_file))

wes <- read.table(wes_file, sep = "\t", header = TRUE) %>%
  dplyr::mutate(CHROM = paste("chr", CHROM, sep = ""))
wgs <- read.table(wgs_file, sep = "\t", header = TRUE) %>%
  dplyr::mutate(CHROM = paste("chr", CHROM, sep = ""))
segdup <- read.delim(segdup_file, sep = "\t", header = TRUE)

combined_data <- ldply(list(WES = wes, WGS = wgs)) 
combined_data$CHROM <- paste("chr", combined_data$CHROM, sep = "")
names(combined_data)[1] <- "Dataset"

# Create GRanges objects ----
wes_grange <- makeGRangesFromDataFrame(
  df = wes,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqnames.field = "CHROM",
  start.field = "POS",
  end.field = "POS",
  starts.in.df.are.0based = FALSE
)
wgs_grange <- makeGRangesFromDataFrame(
  df = wgs,
  keep.extra.columns = TRUE,
  ignore.strand = TRUE,
  seqnames.field = "CHROM",
  start.field = "POS",
  end.field = "POS",
  starts.in.df.are.0based = FALSE
)

# Convert to hg19 ----
hub <- AnnotationHub()
chain <- query(hub, 'hg38ToHg19')[[1]] 

wes_hg19 <- liftOver(wes_grange, chain) %>%
  as.data.frame()
wgs_hg19 <- liftOver(wgs_grange, chain) %>%
  as.data.frame()

################################################################################
segdup$color <- laply(segdup$fracMatch, segdup_color)

chromstart_hg19 = 18656000 / 1e6
chromend_hg19   = 21792000 / 1e6
chromstart_hg38 = 18173233 / 1e6
chromend_hg38 = 21437711 / 1e6

wes_depth <- ggplot(wes_hg19, aes(x = start / 1e6, y = MEAN_DEPTH)) +
  geom_point(size = 1.2) +
  xlim(chromstart_hg19, chromend_hg19) + 
  ylim(0, 300) +
  theme_bw() +
  theme(plot.margin = unit(c(0,0.05,0,0.5), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8),
        title = element_text(size = 10),
        text = element_text(family = "Arial")) +
  ggtitle("WES average depth per site") +
  xlab("Position (Mb)") +
  ylab("Average depth")

wgs_depth <- ggplot(wgs_hg19, aes(x = start / 1e6, y = MEAN_DEPTH)) +
  geom_point(size = 1.2) +
  xlim(chromstart_hg19, chromend_hg19) +
  theme_bw() +
  ylim(0, 300) +
  theme(plot.margin = unit(c(0,0.05,0,0.5), "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x  = element_text(size = 8),
        title = element_text(size = 10),
        text = element_text(family = "Arial")) +
  ggtitle("WGS average depth per site") +
  xlab("Position (Mb)") +
  ylab("Average depth")

# bed
bedgraph <- ggplot(segdup) +
  geom_segment(aes(x     = chromStart / 1e6, 
                   xend  = chromEnd / 1e6, 
                   y     = 100, 
                   yend  = 100,
                   color = factor(color)),
               size = 10,
               alpha = 0.2) +
  scale_color_manual(values = c("1" = "#000000", "2" = "#FFB90F",
                                "3" = "#FF7F00")) +
  ggplot2::xlim(chromstart_hg19, chromend_hg19) +
  ggplot2::ylim(0, 200) +
  theme_bw() +
  theme(plot.margin = unit(c(0,0.05,0,0.5), "cm")) +
  ggtitle("Segmental duplications") +
  xlab("Position (Mb)") +
  ylab("Foo") +
  theme(axis.text.y = element_text(color = "white", size = 8),
        axis.title.y = element_text(color = "white", size = 10),
        axis.ticks.y = element_line(color = "white"),
        axis.text.x  = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        title = element_text(size = 10),
        legend.position = "none",
        text = element_text(family = "Arial"))


# Multiplot
pdf(file = output_file, width = 7.2, height = 3.5)
multiplot(wes_depth, wgs_depth, bedgraph)
dev.off()
