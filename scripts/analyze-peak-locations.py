#!/usr/bin/env python2
import argparse
import gff
import os

# In below three functions:
#   * strand should be +, -, or ., as per GFF spec.
#   * start and end refer to coordinates on chromosome, with start <= end.
#   * upstream and downstream refer to number of nucleotides in respective
#     direction relative to feature, with strandedness taken into account
#     (meaning that on - strand, upstream of feature is farther to right
#     [i.e., larger] in chromosome coordinates).
def calculate_proximal_tss_coords(start, end, strand, upstream, downstream):
  if strand in ('+', '.'):
    left = start - upstream
    right = start + downstream
  elif strand == '-':
    left = end - downstream
    right = end + upstream
  else:
    raise Exception('Unknown strand %s' % strand)
  return (left, right)

def calculate_proximal_tes_coords(start, end, strand, upstream, downstream):
  if strand in ('+', '.'):
    left = end - upstream
    right = end + downstream
  elif strand == '-':
    left = start - downstream
    right = start + upstream
  else:
    raise Exception('Unknown strand %s' % strand)
  return (left, right)

def calculate_inside_gene_coords(
    start, end, strand,
    tss_upstream, tss_downstream,
    tes_upstream, tes_downstream
):
  if strand in ('+', '.'):
    left = start + tss_downstream + 1
    right = end - tes_upstream - 1
  elif strand == '-':
    left = start + tes_upstream + 1
    right = end - tss_downstream - 1
  else:
    raise Exception('Unknown strand %s' % strand)
  return (left, right)

def prepare_reference_sets(transcripts_filename, proximal_promoter_boundaries,
                           proximal_gene_end_boundaries):
  reference_sets = {
    'proximal_tss': {},
    'proximal_tes': {},
    'inside_gene': {},
  }

  proximal_tss_upstream_nts = proximal_promoter_boundaries[0]
  proximal_tss_downstream_nts = proximal_promoter_boundaries[1]
  proximal_tes_upstream_nts = proximal_gene_end_boundaries[0]
  proximal_tes_downstream_nts = proximal_gene_end_boundaries[1]

  transcripts_file = open(transcripts_filename)
  # Remember, in GFF:
  #   * Ranges are indexed starting at 1
  #   * Ranges are inclusive -- i.e., they are defined as [start, end]
  for line in transcripts_file:
    transcript = gff.parse_gff_line(line)
    seqname = transcript['seqname'].lower()
    start = transcript['start']
    end = transcript['end']
    strand = transcript['strand']

    for set_name in ('proximal_tss', 'proximal_tes', 'inside_gene'):
      if seqname not in reference_sets[set_name]:
        reference_sets[set_name][seqname] = []

    reference_sets['proximal_tss'][seqname].append(
      calculate_proximal_tss_coords(start, end, strand,
        proximal_tss_upstream_nts, proximal_tss_downstream_nts
      )
    )
    reference_sets['proximal_tes'][seqname].append(
      calculate_proximal_tes_coords(start, end, strand,
        proximal_tes_upstream_nts, proximal_tes_downstream_nts
      )
    )
    reference_sets['inside_gene'][seqname].append(
      calculate_inside_gene_coords(start, end, strand,
        proximal_tss_upstream_nts, proximal_tss_downstream_nts,
        proximal_tes_upstream_nts, proximal_tes_downstream_nts
      )
    )

  transcripts_file.close()
  return reference_sets

def analyze_peak(peak_line, reference_sets, output_files):
  peak = gff.parse_gff_line(peak_line)
  if peak is None:
    return

  peak_range = (peak['start'], peak['end'])
  seqname = peak['seqname'].lower()

  for range_type in ('proximal_tss', 'proximal_tes', 'inside_gene'):
    if seqname not in reference_sets[range_type]:
      continue
    for location_range in reference_sets[range_type][seqname]:
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
  parser = argparse.ArgumentParser(description='Analyze peak locations.')
  parser.add_argument('peaks_file', help='Path to GFF-formatted peaks file')
  parser.add_argument('transcripts_file', help='Path to transcripts file')
  parser.add_argument('output_dir', help='Path to ouput directory')
  parser.add_argument('--proximal-tss-upstream-nts', type=int, default=500,
    help='Number of nucleotides upstream of TSS that define proximal promoter boundary')
  parser.add_argument('--proximal-tss-downstream-nts', type=int, default=1500,
    help='Number of nucleotides downstream of TSS that define proximal promoter boundary')
  parser.add_argument('--proximal-tes-upstream-nts', type=int, default=1500,
    help='Number of nucleotides upstream of transcription end site that define gene end boundary')
  parser.add_argument('--proximal-tes-downstream-nts', type=int, default=500,
    help='Number of nucleotides downstream of transcription end site that define gene end boundary')
  return parser.parse_args()

def main():
  args = parse_args()
  reference_sets = prepare_reference_sets(args.transcripts_file,
    (args.proximal_tss_upstream_nts, args.proximal_tss_downstream_nts),
    (args.proximal_tes_upstream_nts, args.proximal_tes_downstream_nts),
  )
  analyze_peak_locations(args.peaks_file, reference_sets, args.output_dir)

if __name__ == '__main__':
  main()
