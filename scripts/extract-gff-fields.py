#!/usr/bin/env python2

# Extract fields of interest from GFF file for statistics calculations
# Original author: Mark Bieda
# Rewritten by: Jeff Wintersinger

import argparse
import os

def parse_gff_line(line):
  line = line.strip()
  if line.startswith('#'):
    return None

  tokens = line.split()
  fields = ('seqname', 'source', 'feature', 'start', 'end',
            'score', 'strand', 'frame', 'group')
  parsed = {}

  for token_index in range(len(fields)):
    field_name = fields[token_index]
    parsed[field_name] = tokens[token_index]
  for numeric_field_name in ('start', 'end'):
    parsed[numeric_field_name] = int(parsed[numeric_field_name])
  parsed['score'] = float(parsed['score'])

  return parsed

def write_line(dest_file, token):
  dest_file.write('%s\n' % token)

def extract_fields(gff_path, output_dir):
  output_types = (
    'scores',
    'feature_types',
    'feature_spacings',
    'feature_lengths',
  )
  output_files = {}
  for output_type in output_types:
    path = os.path.join(output_dir, output_type)
    output_files[output_type] = open(path, 'w')
  print output_files

  gff_input = open(gff_path)
  feature_types = set()
  prev_feature_start = None

  for input_line in gff_input:
    parsed_gff = parse_gff_line(input_line)
    if parsed_gff is None:
      continue

    feature_types.add(parsed_gff['feature'])
    feature_length = parsed_gff['end'] - parsed_gff['start'] + 1
    write_line(output_files['feature_lengths'], feature_length)
    write_line(output_files['scores'], parsed_gff['score'])

    current_feature_start = parsed_gff['start']
    if prev_feature_start is not None:
      feature_spacing = current_feature_start - prev_feature_start
      write_line(output_files['feature_spacings'], feature_spacing)
    prev_feature_start = current_feature_start
    
  for feature_type in feature_types:
    write_line(output_files['feature_types'], feature_type)

  gff_input.close()
  for output_type in output_files.keys():
    output_files[output_type].close()

def main():
  parser = argparse.ArgumentParser(description='Extract fields of interest from GFF file for statistics calculations')
  parser.add_argument('gff_path', help='Path to GFF input file')
  parser.add_argument('output_path', help='Path to directory in which to write output')
  args = parser.parse_args()
  extract_fields(args.gff_path, args.output_path)

if __name__ == '__main__':
  main()
