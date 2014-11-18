from collections import defaultdict
from itertools import count
import csv

#Generator of infinite powers of two
def indexer():
  for index in count(0):
    yield pow(2, index)

#Sidecount algorithm, used later to weight adding gene sets
def sidecount(bits):
  for index in count(1):
    bits = (bits - 1) & bits
    if not bits:
      return index

index_generator = indexer()

#data structures
gene_index_lookup = defaultdict(index_generator.next)
index_gene_lookup = dict()
file_hash_lookup = defaultdict(lambda: 0)

#gene_universe = set()

with open("confused_protein_files.txt", 'r') as protein_file:
  reader = csv.reader(protein_file, delimiter='\t')
  reader.next() #skip the header line
  for (gene, filename) in reader:
    #iterate through the lines of the file and load our lookup structures
    index_gene_lookup[gene_index_lookup[gene]] = gene
    file_hash_lookup[filename] += gene_index_lookup[gene]

print "Gene sets loaded."

minimum_set = set()

#So now we have a bijective map of file_name -> binary hash, where every
#nth bit of the hash is the presence or absence of the nth gene. This makes
#it pretty fast to determine whether a file adds new proteins to the min set
#or, better yet, obviates a file we've already added.

universe = index_generator.next() -1 #bithash for the whole gene set
total = 0
best_pass_candidate = None
while total < universe:
  for filename, gene_hash in file_hash_lookup.items():
    if total | gene_hash != total: #if adding it would add new genes
      if sidecount(file_hash_lookup[best_pass_candidate]) <= sidecount(gene_hash):
        best_pass_candidate = filename
  minimum_set.add(best_pass_candidate)
  best_pass_candidate = None
  for f in minimum_set:
    total = total | file_hash_lookup[f]


#well, here's our minimum set, let's spool it out
#in a readable format
def decompose(bithash):
  "decompose integer into a list of powers of 2"
  start = 0
  return_list = list()
  for bit in list(reversed(bin(bithash)))[:-2]: #skip the '0b' at the begnning of the string
    if int(bit):
      return_list.append(pow(2, start))
    start += 1
  return return_list

#print bin(reduce(lambda x, y: y + x, index_gene_lookup.keys()))
for f in minimum_set:
  print f, ":", ','.join([index_gene_lookup[i] for i in decompose(file_hash_lookup[f])])
print len(minimum_set), "files representing", sidecount(universe), "genes."
