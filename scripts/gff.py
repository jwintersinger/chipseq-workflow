def parse_group(group):
  attribs = {}
  for field in group.split(';'):
    if '=' not in field:
      continue
    key, value = field.split('=', 2)
    attribs[key.strip()] = value.strip()
  return attribs

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

  # Make 'start' and 'end' numeric.
  for intergral_field_name in ('start', 'end'):
    parsed[intergral_field_name] = int(parsed[intergral_field_name])

  # Munge score.
  if parsed['score'] == '.':
    parsed['score'] = None
  else:
    parsed['score'] = float(parsed['score'])

  # Parse group.
  parsed['group'] = parse_group(parsed['group'])

  return parsed
