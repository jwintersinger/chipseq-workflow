#!/usr/bin/env Rscript
calculate_score_stats <- function(scores_list_path, scores_stats_output_path) {
  scores <- read.table(scores_list_path)$V1
  min <- min(scores)
  max <- max(scores)
  mean <- mean(scores)
  median <- median(scores)
  p05 <- quantile(scores, 0.05)
  p50 <- quantile(scores, 0.50)
  p90 <- quantile(scores, 0.90)
  p95 <- quantile(scores, 0.95)
  p98 <- quantile(scores, 0.98)
  p99 <- quantile(scores, 0.99)
  numLines <- format(length(scores), big.mark=",")

  sink(scores_stats_output_path)
  cat(
    c(
      paste("Min: ", min),
      paste("Max: ", max),
      paste("Mean: ", mean),
      paste("Median: ", median),
      paste("90th percentile: ", p90),
      paste("95th percentile: ", p95),
      paste("99th percentile: ", p99),
      paste("Number of lines: ", numLines)
    ), sep="\n"
  )
  sink()
}

graph_feature <- function(feature_path, graph_path) {
  feature <- read.table(feature_path)$V1
  pdf(graph_path)
  hist(feature)
  dev.off()
}

main <- function() {
  calculate_score_stats('scores', 'scores_stats')
  graph_feature('scores', 'scores.pdf')
  graph_feature('feature_lengths', 'feature_lengths.pdf')
  graph_feature('feature_spacings', 'feature_spacings.pdf')
}

main()
