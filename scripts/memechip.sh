#!/bin/sh
# Run MEME workflow for ChIP-seq data.
# TODO: switch to using named rather than positional parameters.

peak_regions="$1"
jaspar_path="$2"
reference_genome="$3"
sequence_width="$4"
sequences_to_sample="$5"
max_run_time="$6"
motifs_to_find="$7"
min_motif_width="$8"
max_motif_width="$9"

peak_seqs=peak_sequences.fa
fastaFromBed -fi "$reference_genome" -bed "$peak_regions" -fo "$peak_seqs"

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
