#!/usr/bin/env python2

# Convert BED file output by MACS to GFF file.
# Original author: Mark Bieda
# Rewritten by: Jeff Wintersinger

import argparse

def parse_bed_line(bed_line):
  tokens = bed_line.split()
  # BED specification: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
  attribs = {
    'seq_name': tokens[0],
    'start': int(tokens[1]),
    'end': int(tokens[2]),
    'name': tokens[3],
    'score': float(tokens[4]),
  }
  return attribs

def make_gff_line(bed_attribs):
  # GFF specification: http://genome.ucsc.edu/FAQ/FAQformat.html#format3
  gff_attribs = (
    bed_attribs['seq_name'],  # Chromosome/scaffold name
    'MACS',                   # Program that generated feature
    'bed2gff',                # Type of feature -- probably not a good value
    bed_attribs['start'] + 1, # Start position of feature. BED is 0-based and GFF is 1-based, so add 1.
    bed_attribs['end']    ,   # End position of feature. Note that though GFF is inclusive and BED is exclusive,
                              #   no adjustment is necessary because of modification to start position made above.
    bed_attribs['score'],     # Score
    '.',                      # Strand -- "." means don't know/care
    '.',                      # Reading frame of first base (0-2) -- "." means that this isn't coding exon,
                              #   which may or may not be correct, but given our ignorance as to whether it is,
                              #   this seems the best value.
    bed_attribs['name'],      # Group -- all items with same group linked into single item
  )

  gff_attribs = [str(elem) for elem in gff_attribs]
  gff_line = '\t'.join(gff_attribs)
  return gff_line

def convert(bed_path, gff_path):
  bed_file = open(bed_path)
  gff_file = open(gff_path, 'w')

  for bed_line in bed_file:
    bed_attribs = parse_bed_line(bed_line)
    gff_line = make_gff_line(bed_attribs)
    gff_file.write(gff_line + '\n')

  gff_file.close()
  bed_file.close()

def main():
  parser = argparse.ArgumentParser(description="Convert MACS' BED format to GFF")
  parser.add_argument('bed_file', nargs=1, help='Path to BED file')
  parser.add_argument('gff_file', nargs=1, help='Path to GFF file')
  args = parser.parse_args()

  convert(args.bed_file[0], args.gff_file[0])

if __name__ == '__main__':
  main()
