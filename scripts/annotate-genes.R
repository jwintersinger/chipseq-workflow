# Annotate genes via Bioconductor's ChIPpeakAnno. For each peak, find the gene
# whose TSS is closest, then determine its expression level.
library(ChIPpeakAnno)
library(rtracklayer)

cmdArgs <- commandArgs(trailingOnly = TRUE)
peaksFileNamePrefix <- cmdArgs[1]
ensemblToAffyMapFile <- cmdArgs[2]
expressionLevelsFile <- cmdArgs[3]
cellLine <- cmdArgs[4]
expressionAnalysisThreshold <- as.numeric(cmdArgs[5])
outputDir <- cmdArgs[6]
ensemblMapperScript <- cmdArgs[7]
allPeaksOutputFile <- paste(outputDir, '/allPeaks.gff', sep='')
overlappingPeaksOutputFile <- paste(outputDir, '/overlappingPeaks.gff', sep='')
nonOverlappingPeaksOutputFile <- paste(outputDir, '/nonOverlappingPeaks.gff', sep='')

#peaksName <- '../macs/lin9_peaks'
#ensemblToAffyMap <- '../../data/expression/mapped-names'
#expressionLevels <- '../../data/expression/U133AGNF1B.gcrma.avg.csv'
#cellLine <- 'CD4+_Tcells'
#expressionAnalysisThreshold <- 10
#allPeaksOutput <- '../chippeakanno/all_peaks.gff'
#overlappingPeaksOutput <- '../chippeakanno/overlapping_peaks.gff'
#nonOverlappingPeaksOutput <- '../chippeakanno/non_overlapping_peaks.gff'

data(TSS.human.NCBI36)
peaksBasicName <- paste(peaksFileNamePrefix, '_peaks.gff', sep='')
peaksFullName <- paste(peaksFileNamePrefix, '_peaks.xls', sep='')

peaksBasic <- import(peaksBasicName)
peaksFull <- read.table(peaksFullName, header=T)

# Annotate peaks
annotatedPeaks <- annotatePeakInBatch(peaksBasic, AnnotationData=TSS.human.NCBI36)
# Sort peaks by "peak" column
annotatedPeaks <- annotatedPeaks[with(annotatedPeaks, order(peak)),]
# Add score column.
annotatedPeaks$score <- peaksBasic$score

# Make graphs
pdf(paste(outputDir, '/graphs.pdf', sep=''))
y <- annotatedPeaks$distancetoFeature[!is.na(annotatedPeaks$distancetoFeature) && annotatedPeaks$fromOverlappingOrNearest == "NearestStart"]
y <- y[y %in% -50000:50000]
hist(y, xlab = 'Distance to nearest TSS', main = '', breaks = 100,xlim = c(min(y) - 100, max(y) + 100))

annotatedPeaksDf <- as.data.frame(annotatedPeaks)
pie(table(annotatedPeaksDf [as.character(annotatedPeaksDf$fromOverlappingOrNearest) == "Overlapping" | (as.character(annotatedPeaksDf$fromOverlappingOrNearest) == "NearestStart" & !annotatedPeaksDf$peak %in% annotatedPeaksDf[as.character(annotatedPeaksDf$fromOverlappingOrNearest) == "Overlapping",]$peak),]$insideFeature))
dev.off()

# Add gene expression data.
cmd <- paste(ensemblMapperScript, ensemblToAffyMapFile, expressionLevelsFile, cellLine, sep = ' ')
exprLevel <- system(cmd, intern = TRUE, input = annotatedPeaks$feature)
annotatedPeaks$exprLevel <- as.numeric(exprLevel)

# Summarize expression levels for genes that overlap peaks and genes that do not.
# Eliminate peaks nearest genes that don't have expression data (i.e., when exprLevel == -1).
peaksNearExpressedGenes <- annotatedPeaks[annotatedPeaks$exprLevel > 0,]
# Select only peaks with scores in expressionAnalysisThreshold%.
peaksNearExpressedGenes <- peaksNearExpressedGenes[peaksNearExpressedGenes$score > quantile(peaksNearExpressedGenes$score, prob = 1 - expressionAnalysisThreshold/100),]
# Possible values of $insideFeature:
#   inside:         peak is entirely within gene
#   includeFeature: gene is entirely within peak
#   overlapStart:   peak overlaps start of gene
#   overlapEnd:     peak overlaps end of gene
#   downstream:     peak is downstream of gene
#   upstream:       peak is upstream of gene
peakSelector <- peaksNearExpressedGenes$insideFeature %in% c('inside', 'includeFeature', 'overlapStart', 'overlapEnd')
overlappingPeaks <- peaksNearExpressedGenes[peakSelector,]
nonOverlappingPeaks <- peaksNearExpressedGenes[!peakSelector,]
mean(overlappingPeaks$exprLevel)
mean(nonOverlappingPeaks$exprLevel)
var(overlappingPeaks$exprLevel)
var(nonOverlappingPeaks$exprLevel)

export(annotatedPeaks, allPeaksOutputFile, 'gff3')
export(overlappingPeaks, overlappingPeaksOutputFile, 'gff3')
export(nonOverlappingPeaks, nonOverlappingPeaksOutputFile, 'gff3')
