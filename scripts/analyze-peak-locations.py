#!/usr/bin/env python2
import argparse
import gff
import os

def prepare_reference_sets(transcripts_filename):
  reference_sets = {
    'proximal_tss': [],
    'proximal_tes': [],
    'inside_gene': [],
  }

  proximal_tss_upstream_nts = -1500
  proximal_tss_downstream_nts = 500
  proximal_tes_upstream_nts = -2000
  proximal_tes_downstream_nts = 300

  transcripts_file = open(transcripts_filename)
  # Remember, in GFF:
  #   * Ranges are indexed starting at 1
  #   * Ranges are inclusive -- i.e., they are defined as [start, end]
  for line in transcripts_file:
    transcript = gff.parse_gff_line(line)
    start = transcript['start']
    end = transcript['end']

    reference_sets['proximal_tss'].append((
      start + proximal_tss_upstream_nts,
      start + proximal_tss_downstream_nts
    ))
    reference_sets['proximal_tss'].append((
      end + proximal_tes_upstream_nts,
      end + proximal_tes_downstream_nts
    ))
    reference_sets['inside_gene'].append((
      start + proximal_tss_downstream_nts + 1,
      end + proximal_tes_upstream_nts - 1
    ))

  transcripts_file.close()
  return reference_sets

def analyze_peak(peak_line, reference_sets, output_files):
  peak = gff.parse_gff_line(peak_line)
  peak_range = (peak['start'], peak['end'])

  for range_type in ('proximal_tss', 'proximal_tes', 'inside_gene'):
    for location_range in reference_sets[range_type]:
      if ranges_overlap(location_range, peak_range):
        output_files[range_type].write(peak_line)
        return
  else:
    output_files['distal'].write(peak_line)

def analyze_peak_locations(peaks_filename, reference_sets, output_dir):
  peaks_file = open(peaks_filename)
  output_files = {
    'proximal_tss': 'proximal_tss.gff',
    'proximal_tes': 'proximal_tes.gff',
    'inside_gene': 'inside_gene.gff',
    'distal': 'distal.gff'
  }
  for output_type in output_files.keys():
    path = os.path.join(output_dir, output_files[output_type])
    output_files[output_type] = open(path, 'w')

  for peak_line in peaks_file:
    analyze_peak(peak_line, reference_sets, output_files)
    
  peaks_file.close()
  for output_type in output_files.keys():
    output_files[output_type].close()

def ranges_overlap(a, b):
  if not (len(a) == len(b) == 2):
    raise Exception('Range is not two-tuple')
  x1, y1 = a[0], a[1]
  x2, y2, = b[0], b[1]
  return x1 <= y2 and x2 <= y1
 
def parse_args():
  parser = argparse.ArgumentParser(description='Analyze peak locations')
  parser.add_argument('peaks_file', nargs=1, help='Path to GFF-formatted peaks file')
  parser.add_argument('transcripts_file', nargs=1, help='Path to transcripts file')
  parser.add_argument('output_dir', nargs=1, help='Path to ouput directory')
  return parser.parse_args()

def main():
  args = parse_args()
  reference_sets = prepare_reference_sets(args.transcripts_file[0])
  analyze_peak_locations(args.peaks_file[0], reference_sets, args.output_dir[0])

if __name__ == '__main__':
  main()
