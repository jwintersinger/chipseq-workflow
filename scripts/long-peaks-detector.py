#!/usr/bin/env python2
import argparse
import gff

def detect_long_peaks(input_gff_path, output_gff_path, min_peak_length):
  input_gff = open(input_gff_path)
  output_gff = open(output_gff_path, 'w')

  for line in input_gff:
    parsed_gff = gff.parse_gff_line(line)
    feature_length = parsed_gff['end'] - parsed_gff['start'] + 1
    if feature_length >= min_peak_length:
      output_gff.write(line + "\n")

  output_gff.close()
  output_gff.close()

def main():
  parser = argparse.ArgumentParser(description='Output only peaks with minimum desired length')
  parser.add_argument('input_gff', help='Path to GFF input file')
  parser.add_argument('output_gff', help='Path to GFF output file into which long peaks will be written')
  parser.add_argument('-m', '--min-peak-length', type=int, default=5000,
                      help='Minimum peak length')
  args = parser.parse_args()
  detect_long_peaks(args.input_gff, args.output_gff, args.min_peak_length)

if __name__ == '__main__':
  main()
