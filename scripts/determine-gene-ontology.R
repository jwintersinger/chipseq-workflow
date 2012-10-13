#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

parse_arguments <- function() {
  option_list <- list(
    make_option(c('-o', '--output-dir'),
                dest='output_dir',
                type='character',
                help='Directory in which to store output'),
    make_option(c('-a', '--annotation-db'),
                dest='annotation_db',
                type='character',
                default='org.Hs.eg.db',
                help='Database to use for determing gene ontology information'),
    make_option(c('--p-value-adjuster'),
                dest='p_value_adjuster',
                type='character',
                default='BH',
                help='Method for computing adjusted p-values for multiple testing procedures (see mt.rawp2adjp)'),
    make_option(c('-p', '--max-p-value'),
                dest='max_p_value',
                type='double',
                default=0.05,
                help='Maximum p value for classification to be considered significant'),
    make_option(c('-i', '--gene-ids'),
                dest = 'gene_ids',
                type = 'character',
                help = 'File containing Ensembl gene IDs for which to determine gene ontology')
  )

  option_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(option_parser)

  for(opt in list('output_dir', 'gene_ids')) {
    if(is.null(opts[[opt]])) {
      stop(paste('Missing required argument:', opt))
    }
  }

  return(opts)
}

main <- function() {
  opts <- parse_arguments()

  # Only load ChIPpeakAnno when necessary to prevent unduly extended load
  # times when unneeded.
  suppressPackageStartupMessages(library(ChIPpeakAnno))

  ensembl_gene_ids <- readLines(opts$gene_ids)
  gene_ontology <- getEnrichedGO(
    ensembl_gene_ids,
    orgAnn = opts$annotation_db,
    feature_id_type = 'ensembl_gene_id',
    maxP = opts$max_p_value,
    multiAdj = TRUE,
    multiAdjMethod = opts$p_value_adjuster
  )

  output_filenames = list(
    bp = 'biological_processes',
    mf = 'molecular_functions',
    cc = 'cellular_components'
  )

  for(go_class in c('bp', 'mf', 'cc')) {
    go_data <- gene_ontology[[go_class]]
    # Sort by count.InDataset column.
    go_data <- go_data[with(go_data, order(count.InDataset, decreasing = TRUE)),]
    go_df <- data.frame(
      'entrez_id' = go_data$EntrezID,
      'go_id' = go_data$go.id,
      'go_term' = go_data$go.term,
      'count_in_dataset' = go_data$count.InDataset,
      'count_in_genome' = go_data$count.InGenome,
      'p_value' = go_data$BH.adjusted.p.value
    )
    output_filename = paste(opts$output_dir, '/', output_filenames[[go_class]], '.csv', sep='')
    write.table(go_df, output_filename, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
  }
}

main()
