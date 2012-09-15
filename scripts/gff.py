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
