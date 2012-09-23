# Annotate genes via Bioconductor's ChIPpeakAnno. For each peak, find the gene
# whose TSS is closest, then determine its expression level.
suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(rtracklayer))

# TODO: rewrite to use optparse package in CRAN.
cmdArgs <- commandArgs(trailingOnly = TRUE)
peaksBasicName <- paste(cmdArgs[1])
peaksFullName <- cmdArgs[2]
ensemblToAffyMapFile <- cmdArgs[3]
expressionLevelsFile <- cmdArgs[4]
cellLine <- cmdArgs[5]
expressionAnalysisThreshold <- as.numeric(cmdArgs[6])
outputDir <- cmdArgs[7]
ensemblMapperScript <- cmdArgs[8]
annotationDb <- cmdArgs[9]
referenceAssembly <- cmdArgs[10]

allPeaksOutputFile <- paste(outputDir, '/allPeaks.gff', sep='')
overlappingPeaksOutputFile <- paste(outputDir, '/overlappingPeaks.gff', sep='')
nonOverlappingPeaksOutputFile <- paste(outputDir, '/nonOverlappingPeaks.gff', sep='')

peaksBasic <- import(peaksBasicName)
peaksFull <- read.table(peaksFullName, header=T)

# Annotate peaks
data(list=c(referenceAssembly))
annotatedPeaks <- annotatePeakInBatch(peaksBasic, AnnotationData=get(referenceAssembly))
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

# Add HUGO gene symbols.
annotatedPeaks <- addGeneIDs(annotatedPeaks, annotationDb, c('symbol'))

# Add gene ontology information.
enrichedGO <- getEnrichedGO(annotatedPeaks, orgAnn=annotationDb)
write.csv(enrichedGO$bp, paste(outputDir, '/biological_processes.csv', sep=''))
write.csv(enrichedGO$cc, paste(outputDir, '/cellcular_components.csv', sep=''))
write.csv(enrichedGO$mf, paste(outputDir, '/molecular_functions.csv', sep=''))

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
