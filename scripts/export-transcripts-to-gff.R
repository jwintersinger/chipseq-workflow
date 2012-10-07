#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))

parse_arguments <- function() {
  option_list <- list(
    # Other possible value: TSS.human.NCBI36
    make_option(c('-r', '--reference-assembly'),
                dest='reference_assembly',
                type='character',
                default='TSS.human.GRCh37',
                help='Reference assembly to use for peak annotation')
  )

  option_parser <- OptionParser(option_list = option_list)
  args <- parse_args(option_parser, positional_arguments = TRUE)

  return(args)
}

massage_chromsome_name <- function(x) {
  if(x %in% c(1:23, 'X', 'Y')) {
    return(paste('chr', x, sep=''))
  } else {
    # I don't know why this must be converted to a string -- otherwise, it
    # returns a nonsensical numerical value.
    return(toString(x))
  }
}

main <- function() {
  args <- parse_arguments()
  if(length(args$args) != 1) {
    stop('Wrong number of required positional arguments')
  } else {
    output_filename <- args$args[1]
  }

  tss_list_name <- args$options$reference_assembly
  # Delay loading ChIPpeakAnno until absolutely necessary given its long load
  # time. Otherwise, invocations of the script that fail take a long time to
  # error out.
  suppressPackageStartupMessages(library(ChIPpeakAnno))
  data(list=c(tss_list_name))
  tss_list <- get(tss_list_name)
  tss_df <- as.data.frame(tss_list)

  tss_df$space <- sapply(tss_df$space, massage_chromsome_name)
  gff_df <- data.frame(tss_df$space, tss_list_name, tss_df$names,
                       tss_df$start, tss_df$end, '.', tss_df$strand)
  write.table(gff_df, output_filename, sep='\t', quote=FALSE, col.names=FALSE, row.names=FALSE)
}

main()
