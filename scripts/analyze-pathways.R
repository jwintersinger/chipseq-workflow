#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

parse_arguments <- function() {
  option_list <- list(
    make_option(c('-o', '--output-file'),
                dest='output_filename',
                type='character',
                help='File in which to store output'),
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
    make_option(c('-m', '--min-term-count'),
                dest='min_term_count',
                type='integer',
                default=10,
                help='Minimum count for pathway term in genome to be included'),
    make_option(c('-i', '--gene-ids'),
                dest = 'gene_ids',
                type = 'character',
                help = 'File containing Ensembl gene IDs for which to determine gene ontology')
  )

  option_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(option_parser)

  for(opt in list('output_filename', 'gene_ids')) {
    if(is.null(opts[[opt]])) {
      stop(paste('Missing required argument:', opt))
    }
  }

  return(opts)
}

main <- function() {
  opts <- parse_arguments()

  # Only load ChIPpeakAnno and reactome.db when necessary to prevent unduly
  # extended load times when unneeded.
  suppressPackageStartupMessages(library(ChIPpeakAnno))
  suppressPackageStartupMessages(library(reactome.db))

  ensembl_gene_ids <- readLines(opts$gene_ids)
  # Use try() so that if error relating to "No enriched pathway can be found"
  # is produced, the script does not exit with a non-zero error code.
  tryCatch({
    pathways <- getEnrichedPATH(
      ensembl_gene_ids,
      orgAnn = opts$annotation_db,
      pathAnn = 'reactome.db',
      feature_id_type = 'ensembl_gene_id',
      maxP = opts$max_p_value,
      minPATHterm = opts$min_term_count,
      multiAdjMethod = opts$p_value_adjuster
    )},
    error = function(err) {
      # Output empty file to indicate no pathways found.
      cat('', file=opts$output_filename)
      print(err)
      quit()
    }
  )

  # Sort by count.InDataset column.
  pathways <- pathways[with(pathways, order(count.InDataset, decreasing = TRUE)),]
  pathways_df <- data.frame(
    'entrez_id' = pathways$EntrezID,
    'kegg_path_id' = pathways$path.id,
    'count_in_dataset' = pathways$count.InDataset,
    'count_in_genome' = pathways$count.InGenome,
    'p_value' = pathways$BH.adjusted.p.value
  )
  write.table(pathways_df, opts$output_filename, sep='\t', quote=FALSE, col.names=TRUE, row.names=FALSE)
}

main()
