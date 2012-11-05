#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

parse_arguments <- function() {
  option_list <- list(
    make_option(c('-b', '--basic-peaks'),
                dest='basic_peaks',
                type='character',
                help='Path to GFF file containing basic peaks data'),
    make_option(c('-f', '--extended-peaks'),
                dest='extended_peaks',
                type='character',
                help='Path to CSV file containing extended peaks data'),
    make_option(c('-o', '--output'),
                dest='output_file',
                type='character',
                help='File in which to store annotated peaks'),
    make_option(c('-a', '--annotation-db'),
                dest='annotation_db',
                type='character',
                default='org.Hs.eg.db',
                help='Database to use for determing gene symbols and ontology information'),
    make_option(c('-r', '--reference-assembly'),
                dest='reference_assembly',
                type='character',
                default='TSS.human.GRCh37',
                help='Reference assembly to use for peak annotation')
  )

  option_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(option_parser)

  for(opt in list('basic_peaks', 'extended_peaks', 'output_file')) {
    if(is.null(opts[[opt]])) {
      stop(paste('Missing required argument:', opt))
    }
  }

  return(opts)
}

massage_chromosome_name <- function(chrom) {
  return(ifelse(
    chrom %in% c(1:23, 'X', 'Y'),
    paste('chr', chrom, sep=''),
    # I don't know why this must be converted to a string -- otherwise, it
    # returns a nonsensical numerical value.
    toString(chrom)
  ));
}

annotate_peaks <- function(
  basic_peaks_filename,
  extended_peaks_filename,
  reference_assembly,
  annotation_db,
  feature_location # 'TSS' or 'geneEnd'
) {
  peaks_basic <- import.gff(basic_peaks_filename)
  peaks_extended <- read.table(extended_peaks_filename, header=T)

  # Annotate peaks
  data(list=c(reference_assembly))
  annotated_peaks <- annotatePeakInBatch(
    peaks_basic,
    PeakLocForDistance = 'middle',
    FeatureLocForDistance = 'TSS',
    AnnotationData = get(reference_assembly)
  )
  # Sort peaks by "peak" column
  annotated_peaks <- annotated_peaks[with(annotated_peaks, order(peak)),]
  # Add score column.
  annotated_peaks$score <- peaks_basic$score
  # Add HUGO gene symbols.
  annotated_peaks <- addGeneIDs(annotated_peaks, annotation_db, c('symbol'))

  # Munge data
  annotated_peaks <- as.data.frame(annotated_peaks)
  # Drop unneeded columns.
  to_drop <- c('names', 'peak', 'fromOverlappingOrNearest')
  annotated_peaks <- annotated_peaks[,!(names(annotated_peaks) %in% to_drop)]
  # Rename poorly-named columns.
  names(annotated_peaks)[names(annotated_peaks) == 'feature'] <- 'gene'
  names(annotated_peaks)[names(annotated_peaks) == 'start_position'] <- 'gene_start'
  names(annotated_peaks)[names(annotated_peaks) == 'end_position'] <- 'gene_end'
  # Prepend "chr" to chromosome names.
  annotated_peaks$space <- massage_chromosome_name(annotated_peaks$space)

  return(annotated_peaks)
}

main <- function() {
  opts <- parse_arguments()

  # Only load these libraries when necessary to prevent unduly extended load
  # times when unneeded.
  suppressPackageStartupMessages(library(ChIPpeakAnno))
  suppressPackageStartupMessages(library(rtracklayer))

  annotated_peaks <- annotate_peaks(
    opts$basic_peaks,
    opts$extended_peaks,
    opts$reference_assembly,
    opts$annotation_db,
    'TSS'
  )

  export.gff3(annotated_peaks, opts$output_file)
}

main()
