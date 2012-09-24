# Annotate genes via Bioconductor's ChIPpeakAnno. For each peak, find the gene
# whose TSS is closest, then determine its expression level.

suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(rtracklayer))
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
    make_option(c('-e', '--ensembl-to-affy-map'),
                dest='ensembl_to_affy_map',
                type='character',
                help='File mapping Ensembl to Affy IDs'),
    make_option(c('-m', '--ensembl-mapper-script'),
                dest='ensembl_mapper_script',
                type='character',
                help='Script to map Ensembl to Affy IDs, using file specified via --ensembl-to-affy-map'),
    make_option(c('-l', '--expression-levels'),
                dest='expression_levels',
                type='character',
                help='File containing gene expression levels data'),
    make_option(c('-c', '--cell-line'),
                dest='cell_line',
                type='character',
                help='Name of cell line in expression levels file to use for gene expression data'),
    make_option(c('-t', '--expression-threshold'),
                dest='expression_threshold',
                type='double',
                default=10,
                help='Threshold to use when selecting percentile of peaks for overlapping/nonoverlapping'),
    make_option(c('-o', '--output-dir'),
                dest='output_dir',
                type='character',
                help='Directory in which to store output'),
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

  for(opt in list('basic_peaks', 'extended_peaks', 'ensembl_to_affy_map',
                  'ensembl_mapper_script', 'expression_levels', 'cell_line',
                  'output_dir')) {
    if(is.null(opts[[opt]])) {
      stop(paste('Missing required argument:', opt))
    }
  }

  return(opts)
}

annotate_peaks <- function(basic_peaks_filename, extended_peaks_filename, reference_assembly, annotation_db) {
  peaks_basic <- import(basic_peaks_filename)
  peaks_extended <- read.table(extended_peaks_filename, header=T)

  # Annotate peaks
  data(list=c(reference_assembly))
  annotated_peaks <- annotatePeakInBatch(peaks_basic, AnnotationData=get(reference_assembly))
  # Sort peaks by "peak" column
  annotated_peaks <- annotated_peaks[with(annotated_peaks, order(peak)),]
  # Add score column.
  annotated_peaks$score <- peaks_basic$score
  # Add HUGO gene symbols.
  print(annotation_db)
  annotated_peaks <- addGeneIDs(annotated_peaks, annotation_db, c('symbol'))

  return(annotated_peaks)
}

generate_histograms <- function(annotated_peaks, output_dir) {
  pdf(paste(output_dir, '/graphs.pdf', sep=''))

  y <- annotated_peaks$distancetoFeature[!is.na(annotated_peaks$distancetoFeature) && annotated_peaks$fromOverlappingOrNearest == "NearestStart"]
  y <- y[y %in% -50000:50000]
  hist(y, xlab = 'Distance to nearest TSS', main = '', breaks = 100,xlim = c(min(y) - 100, max(y) + 100))
  
  annotated_peaks_df <- as.data.frame(annotated_peaks)
  pie(table(annotated_peaks_df [as.character(annotated_peaks_df$fromOverlappingOrNearest) == "Overlapping" | (as.character(annotated_peaks_df$fromOverlappingOrNearest) == "NearestStart" & !annotated_peaks_df$peak %in% annotated_peaks_df[as.character(annotated_peaks_df$fromOverlappingOrNearest) == "Overlapping",]$peak),]$insideFeature))

  dev.off()
}

write_gene_ontology_info <- function(annotated_peaks, output_dir, annotation_db) {
  # Add gene ontology information.
  print(annotation_db)
  enriched_go <- getEnrichedGO(annotated_peaks, orgAnn=annotation_db)
  write.csv(enriched_go$bp, paste(output_dir, '/biological_processes.csv', sep=''))
  write.csv(enriched_go$cc, paste(output_dir, '/cellcular_components.csv', sep=''))
  write.csv(enriched_go$mf, paste(output_dir, '/molecular_functions.csv', sep=''))
}


add_gene_expression_data <- function(ensembl_mapper_script, ensembl_to_affy_map,
                                     expression_levels_file, cell_line, annotated_peaks) {
  # Add gene expression data.
  cmd <- paste(ensembl_mapper_script, ensembl_to_affy_map,
               expression_levels_file, cell_line, sep = ' ')
  expr_level <- system(cmd, intern = TRUE, input = annotated_peaks$feature)
  annotated_peaks$exprLevel <- as.numeric(expr_level)
  return(annotated_peaks)
}

perform_location_analysis <- function(annotated_peaks, expression_threshold) {
  # Summarize expression levels for genes that overlap peaks and genes that do not.
  # Eliminate peaks nearest genes that don't have expression data (i.e., when exprLevel == -1).
  peaks_near_expressed_genes <- annotated_peaks[annotated_peaks$exprLevel > 0,]
  # Select only peaks with scores in expression_threshold%.
  peaks_near_expressed_genes <- peaks_near_expressed_genes[peaks_near_expressed_genes$score > quantile(peaks_near_expressed_genes$score, prob = 1 - expression_threshold/100),]

  # Possible values of $insideFeature:
  #   inside:         peak is entirely within gene
  #   includeFeature: gene is entirely within peak
  #   overlapStart:   peak overlaps start of gene
  #   overlapEnd:     peak overlaps end of gene
  #   downstream:     peak is downstream of gene
  #   upstream:       peak is upstream of gene
  peak_selector <- peaks_near_expressed_genes$insideFeature %in% c('inside', 'includeFeature', 'overlapStart', 'overlapEnd')
  overlapping_peaks <- peaks_near_expressed_genes[peak_selector,]
  non_overlapping_peaks <- peaks_near_expressed_genes[!peak_selector,]

  # Currently unused.
  mean(overlapping_peaks$exprLevel)
  mean(non_overlapping_peaks$exprLevel)
  var(overlapping_peaks$exprLevel)
  var(non_overlapping_peaks$exprLevel)

  return(list(overlapping=overlapping_peaks, non_overlapping=non_overlapping_peaks))
}

main <- function() {
  opts <- parse_arguments()


  annotated_peaks <- annotate_peaks(opts$basic_peaks, opts$extended_peaks,
                                    opts$reference_assembly, opts$annotation_db)
  generate_histograms(annotated_peaks, opts$output_dir)
  write_gene_ontology_info(annotated_peaks, opts$output_dir, opts$annotation_db)
  annotated_peaks <- add_gene_expression_data(opts$ensembl_mapper_script, opts$ensembl_to_affy_map,
                           opts$expression_levels, opts$cell_line, annotated_peaks)
  segregated_peaks <- perform_location_analysis(annotated_peaks, opts$expression_threshold)


  all_peaks_output <- paste(opts$output_dir, '/allPeaks.gff', sep='')
  overlapping_peaks_output <- paste(opts$output_dir, '/overlappingPeaks.gff', sep='')
  non_overlapping_peaks_output <- paste(opts$output_dir, '/nonOverlappingPeaks.gff', sep='')

  export_format <- 'gff3'
  export(annotated_peaks, all_peaks_output, export_format)
  export(segregated_peaks$overlapping, overlapping_peaks_output, export_format)
  export(segregated_peaks$non_overlapping, non_overlapping_peaks_output, export_format)
}

main()
