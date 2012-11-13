#!/bin/sh
# Run MEME workflow for ChIP-seq data.

while getopts "p:j:r:w:s:t:f:m:n:o:" opt; do
  case $opt in
    p)
      peak_regions=$OPTARG
      ;;
    j)
      jaspar_path=$OPTARG
      ;;
    r)
      reference_genome=$OPTARG
      ;;
    w)
      sequence_width=$OPTARG
      ;;
    s)
      sequences_to_sample=$OPTARG
      ;;
    t)
      max_run_time=$OPTARG
      ;;
    f)
      motifs_to_find=$OPTARG
      ;;
    m)
      min_motif_width=$OPTARG
      ;;
    n)
      max_motif_width=$OPTARG
      ;;
    o)
      output_dir=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument" >&2
      exit 1
      ;;
  esac
done

for opt in \
  peak_regions jaspar_path reference_genome sequence_width \
  sequences_to_sample max_run_time motifs_to_find \
  min_motif_width max_motif_width; do
  eval 'if [ -z $'${opt}' ]; then echo "Missing parameter: '${opt}'" >&2; exit 1; fi'
done

# Determine FASTA-formatted sequences specified by interval file.
peak_seqs=peak_sequences.fa
fastaFromBed -fi "$reference_genome" -bed "$peak_regions" -fo "$peak_seqs"

# Run MEME suite for ChIP data.
fasta-center -len "$sequence_width" < "$peak_seqs" > ./seqs-centered
fasta-dinucleotide-shuffle -f ./seqs-centered -t -dinuc > ./seqs-shuffled
cat ./seqs-centered ./seqs-shuffled > ./seqs-centered_w_bg
fasta-subsample ./seqs-centered "$sequences_to_sample" -rest ./seqs-discarded > ./seqs-sampled
meme ./seqs-sampled -oc meme_out -dna -mod zoops -nmotifs "$motifs_to_find" -minw "$min_motif_width" -maxw "$max_motif_width" -time "$max_run_time" -revcomp -nostatus
tomtom -verbosity 1 -oc meme_tomtom_out -min-overlap 5 -dist pearson -evalue -thresh 0.1 -no-ssc meme_out/meme.txt "$jaspar_path"
dreme -v 1 -oc dreme_out -p ./seqs-centered -n ./seqs-shuffled -png
tomtom -verbosity 1 -oc dreme_tomtom_out -min-overlap 5 -dist pearson -evalue -thresh 0.1 -no-ssc dreme_out/dreme.txt "$jaspar_path"
ame --verbose 1 --oc ame_out --fix-partition 820 --bgformat 0 ./seqs-centered_w_bg "$jaspar_path"
centrimo -verbosity 1 -oc centrimo_out "$peak_seqs" meme_out/meme.txt dreme_out/dreme.txt
spamo -verbosity 1 -oc meme_spamo_out "$peak_seqs" meme_out/meme.txt meme_out/meme.txt dreme_out/dreme.txt
spamo -verbosity 1 -oc dreme_spamo_out "$peak_seqs" dreme_out/dreme.txt meme_out/meme.txt dreme_out/dreme.txt

# Move output files.
mkdir --parents $output_dir
for dir in {ame,centrimo,dreme,dreme_spamo,dreme_tomtom,meme,meme_spamo,meme_tomtom}_out; do
  mv $dir/* $output_dir/
done

# Echo output.
echo \
  '<!DOCTYPE html>' \
  '<head><title>MEME Results</title></head>' \
  '<body><h1>MEME Results</h1><ul>' \
  '<li>AME: <a href="ame.html">HTML</a>, <a href="ame.txt">text</a></li>' \
  '<li>CentriMo: <a href="centrimo.html">HTML</a>, <a href="centrimo.txt">text</a>, <a href="site_counts.txt">site counts</a></li>' \
  '<li>DREME: <a href="dreme.html">HTML</a>, <a href="dreme.txt">text</a>, <a href="dreme.xml">XML</a></li>' \
  '<li>DREME SpaMo: <a href="spamo.html">HTML</a>, <a href="spamo.xml">XML</a></li>' \
  '<li>DREME TOMTOM: <a href="tomtom.html">HTML</a>, <a href="tomtom.txt">text</a>, <a href="tomtom.xml">XML</a></li>' \
  '<li>MEME: <a href="meme.html">HTML</a>, <a href="meme.txt">text</a>, <a href="meme.xml">XML</a></li>' \
  '<li>MEME SpaMo: <a href="spamo.html">HTML</a>, <a href="spamo.txt">XML</a></li>' \
  '<li>MEME TOMTOM: <a href="tomtom.html">HTML</a>, <a href="tomtom.txt">text</a>, <a href="tomtom.xml">XML</a></li>' \
  '</ul></body></html>' \
  > index.html
