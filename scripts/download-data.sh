#!/bin/sh
# Download all data sets necessary to run workflow.

# Download the human reference genome and concatenate into a single Fasta file.
# Needed to pull out peak sequences for motif finding.

function download_ncbi36() {
  filename=chromFa.zip
  http://128.114.119.163/goldenPath/hg18/bigZips/$filename
  unzip $filename
  rm $filename
  cat chr?.fa chr??.fa > combined.fa
  rm chr*.fa
}

function download_grch37() {
  filename=chromFa.tar.gz
  wget http://128.114.119.163/goldenPath/hg19/bigZips/$filename
  tar xvzf $filename
  rm $filename
  cat chr?.fa chr??.fa > combined.fa
  rm chr*.fa
}

OLD_WD="$(pwd)"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p $SCRIPT_DIR/../data
cd $SCRIPT_DIR/../data

mkdir -p reference-genome/GRCh37 hon sada

# Reference genome
cd reference-genome/GRCh37
download_grch37
cd ..

# Hon's histone methylation data
# Note: this is commented out because the current pipeline is focused on TF binding, not histone methylation.
#cd ../hon
# H3K9me3 data
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM721nnn/GSM721135/GSM721135_HCC1954.merged.H3K9me3.nodup.bam
# Control data
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM721nnn/GSM721139/GSM721139_HCC1954.merged.input.nodup.bam

# Sada's LIN9 data
cd ../sada
# LIN9 data
wget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM665nnn/GSM665907/GSM665907_Lin9_1.bam
# Control data
wget ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples/GSM665nnn/GSM665905/GSM665905_Input_1.bam

cd "$OLD_WD"
