# source("http://bioconductor.org/biocLite.R")
# biocLite("Gviz")

library(Gviz)
library(biomaRt)
library(GenomicFeatures)
library(BSgenome.Hsapiens.NCBI.GRCh37) # Custom BCBI genome sequence annotations
library(rtracklayer)
options(ucscChromosomeNames=FALSE) # To prevent using "chr" prefix in chromosome
bam_file <- "E:/Exome_sequencing_R186/bam/TW122.42421.bam"

genome <- "hg19"
chr <- "22"
pos <- 19455430


# For gene model.
gene_gtf_file <- "C:/Lab/data/genome_wide_gene_list/refseq_genes.hg19.gtf"
gene_features <- GenomicFeatures::makeTxDbFromGFF(file       = gene_gtf_file,
                                                  format     = "gtf",
                                                  organism   = "Homo sapiens",
                                                  dataSource = "RefGene")

# Ideogram ----
ideogram_file <- "C:/Lab/data/genome_wide_gene_list/ucsc.hg19.ideogram.txt"
ideogram_data <- read.table(ideogram_file, sep = "\t", header = TRUE)

# Plot tracks ----
afrom <- pos - 50
ato   <- pos + 50

ideoTrack <- Gviz::IdeogramTrack(bands = ideogram_data, genome = genome, 
                                 chromosome = chr)
axisTrack <- Gviz::GenomeAxisTrack()
sTrack <- SequenceTrack(Hsapiens)

alTrack <- Gviz::AlignmentsTrack(range      = bam_file,
                                 start      = afrom,
                                 end        = ato,
                                 chromosome = chr,
                                 isPaired   = TRUE,
                                 stacking = "dense")

gTrack <- Gviz::GeneRegionTrack(range = gene_features,
                                chromosome = chr,
                                start = afrom,
                                end = ato,
                                stacking = "dense",
                                geneSymbols = TRUE,
                                size = 0.2)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# gTrack <- BiomartGeneRegionTrack(genome      = genome,
#                               chromosome  = chr,
#                               start       = afrom, 
#                               end         = ato, 
#                               biomart     = mart,
#                               stacking    = "full",
#                               geneSymbols = TRUE)
# Plot ----
plotTracks(c(ideoTrack, axisTrack, sTrack, alTrack, gTrack),
           from = afrom, to = ato, chromosome = chr)

