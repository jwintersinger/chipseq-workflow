#!/usr/bin/env python2
import argparse
import gff

def parse_args():
  parser = argparse.ArgumentParser(description='Extract gene names from annotated GFF file.')
  parser.add_argument('gff_file', help='Path to GFF file')
  return parser.parse_args()

def main():
  args = parse_args()
  with open(args.gff_file) as gff_file:
    already_printed = []
    for line in gff_file:
      parsed = gff.parse_gff_line(line)
      if parsed is None:
        continue

      gene_name = parsed['group']['gene']
      if gene_name in already_printed:
        continue

      print gene_name
      already_printed.append(gene_name)

if __name__ == '__main__':
  main()
