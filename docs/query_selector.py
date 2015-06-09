
FATools uses YAML-formated query as its main query system.

Example for samples:

samples:
  S1:
    - { batch: MYBATCH, codes: [ SAMPLE-001, SAMPLE-002 ] }
    - { query: Indonesia[country] }
    - { sample_ids: [ 102, 103, 104, 105] }
  S2:
    - { query: !Indonesia[country] }


Keys:
batch
codes
query
sample_ids

