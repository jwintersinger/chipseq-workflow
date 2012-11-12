#!/bin/bash
while getopts "p:" opt; do
  case $opt in
    p)
      gff_path=$OPTARG
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

if [[ !(-f "$gff_path") || (-z $gff_path)]]; then
  echo "GFF file does not exist: $gff_path" >&2
  exit 1
else
  echo weiner
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Output GFF fields of interest to current directory.
$SCRIPT_DIR/extract-gff-fields.py $gff_path .
$SCRIPT_DIR/generate-peak-stats-and-plots.R 2>/dev/null
