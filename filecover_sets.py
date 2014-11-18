from collections import defaultdict
import csv
#the basic approach.

minimum_set = defaultdict(set)

with open('confused_protein_files.txt', 'r') as protein_file:
  rdr = csv.reader(protein_file, delimiter='\t')
  rdr.next()
  for (gene, filename) in rdr:
    minimum_set[filename].add(gene)

print "Gene sets loaded."

universe = set()
for g in minimum_set.values():
  universe = universe.union(g)

for filename, gene_set in minimum_set.copy().items():
  total_set = set()
  for f, s in minimum_set.items():
    if filename not in f:
      total_set = total_set.union(s)
  if len(total_set) == len(universe):
    if filename in minimum_set:
      del minimum_set[filename]

for filename, gene_set in minimum_set.items():
  print filename, ":", ",".join(gene_set)

print len(minimum_set), "files representing", len(total_set), "genes."
