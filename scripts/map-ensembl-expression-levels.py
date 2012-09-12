#!/usr/bin/env python2

# Map Ensembl gene IDs to corresponding expression levels.
#
# Command line: map-ensembl-expression-levels.py [map file of Ensembl gene IDs to Affymetrix IDs] [gene expression level file]
#               (See Kepler workflow annotations for details on file format.)
# STDIN: list of Ensembl gene IDs, with one per line
# STDOUT: Suppose that line n of STDIN is "GENE".
#   If GENE could be mapped to an Affymetrix U133AGNF1B ID and expression data exists for it:
#     Line n of output will be decimal value representing GENE's expression level
#   Otherwise:
#     Line n of output will be -1

import csv
import os
import sys

def create_ensembl_to_affy_map(map_filename):
  id_map = {}
  with open(map_filename) as mapped_names:
    for line in mapped_names:
      mapped_name = line.split()
      if len(mapped_name) != 2:
        continue
      affy_id, ensembl_id = mapped_name
      id_map[ensembl_id] = affy_id
  return id_map

def map_ensembl_to_affy(ensembl_ids, map_filename):
  ensembl_to_affy_map = create_ensembl_to_affy_map(map_filename)
  affy_ids = []
  for ensembl_id in ensembl_ids:
    ensembl_id = ensembl_id.strip()
    if ensembl_id in ensembl_to_affy_map:
      affy_ids.append(ensembl_to_affy_map[ensembl_id])
    else:
      affy_ids.append('')
  return affy_ids

def fetch_expression_data(expression_data_filename, cell_line):
  cell_line_data = {}
  with open(expression_data_filename) as expression_data:
    reader = csv.DictReader(expression_data)
    for gene in reader:
      affy_gene_id = gene['']
      level = gene[cell_line]
      cell_line_data[affy_gene_id] = level
  return cell_line_data

def main():
  if len(sys.argv) != 4:
    sys.exit('Usage: %s [map file of Ensembl gene IDs to Affymetrix IDs] [gene expression level file] [cell line name]' % sys.argv[0])
  else:
    map_filename, expression_data_filename, cell_line = sys.argv[1:4]

  ensembl_ids = [ensembl_id.strip() for ensembl_id in sys.stdin.readlines()]
  affy_ids = map_ensembl_to_affy(ensembl_ids, map_filename)

  expr_data = fetch_expression_data(expression_data_filename, cell_line)
  for affy_id in affy_ids:
    if affy_id in expr_data:
      print expr_data[affy_id]
    else:
      print -1

if __name__ == '__main__':
  main()
