#!/usr/bin/env python3

'''
Author: Jimmy O'Donnell <jodonnellbio@gmail.com>

source: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc132

Downloading a set of DNA sequences from NCBI's nucleotide database (also known as GenBank or nt) requires two steps:

1. Get a list of sequences that match a query.

2. Download each of the sequences in that list.

For the first step, we use Bio.Entrez.esearch(); for the second, Bio.Entrez.efetch()

Looping over sequences for the download portion would put undo strain on NCBI's servers; therefore it is recommended to get the list of sequences, leave the list on their server, and reference that to download sequences.
To do this, call Bio.Entrez.esearch() as normal, but with the additional argument of usehistory="y"

TODO:
  - query file:
    - columns: "field","include","value"
    - example: "organism", "include", "embiotocidae"
    - example: "organism", "exclude", "ditrema"
    - example: "organism", "exclude", "neoditrema"
'''

# you will need these packages:
from datetime import datetime
from Bio import Entrez
import sys
import argparse

parser = argparse.ArgumentParser(description = 
    'Download sequences from NCBI database nt ("GenBank")')

parser.add_argument('-e', '--email',
  help = 'An email address at which NCBI can contact you in case of overuse.', 
  required = False)

parser.add_argument('-q', '--query',
  help = 'String of search terms, e.g. "Homo sapiens"[Organism]', 
  required = False)

parser.add_argument('-f', '--format',
  choices=['gb', 'fasta'], default = "gb", 
  help = 'Output format (genbank or fasta). Use gb for ecoprimers.', 
  required = False)

parser.add_argument('-o', '--output',
  help = 'Base name for output without extension. ex: oophaga_16s ', 
  required = False)

args = vars(parser.parse_args())

# args will be a dictionary containing the arguments:

if args['email'] is None:
    # you MUST enter an email address:
    print('You must enter an email address. NCBI will contact you in case of overuse.')
    Entrez.email = input("email: ")
else:
    Entrez.email = args['email']

if args['query'] is None:
    print('Enter a search query, for example:' + '\n' + 
    '"Poecilia wingei"[Organism] AND gene_in_mitochondrion[PROP] AND 100:1000[SLEN]')
    my_query = input()
    if len(my_query) < 1:
        print("No query entered. For help, use argument -h.")
        sys.exit()
else:
    my_query = args['query']

# Example queries:
# my_query = '"Embiotocidae"[organism] AND gene_in_mitochondrion[PROP] AND ("cytochrome b") NOT ("COI")'
# my_query = '"Poecilia wingei"[Organism]'
# my_query = '"Boraras brigittae"[Organism]'
# my_query = '"Oophaga pumilio"[Organism] AND ("cytochrome"[ALL] NOT "cytochrome b"[ALL])'

outfile_format = args['format']

start_time = datetime.now()
start_time_fmt = start_time.strftime('%y%m%d-%H%M%S')

if args['output'] is None:
    outfile_base  = 'gbdl' + '_' + start_time_fmt
else:
    outfile_base = args['output']

outfile_path  = outfile_base + '.' + outfile_format
metadata_path = outfile_base + '.md'

# database can be either 'nucleotide' or others, 
# this should probably not be modified by the user of this script
ncbi_db = "nucleotide"

search_handle = Entrez.esearch(db = ncbi_db, usehistory = "y", retmax=200, 
    term = my_query)
search_results = Entrez.read(search_handle)
search_handle.close()


# When you get the XML output back, it will still include the usual search results:
gi_list = search_results["IdList"]
count = int(search_results["Count"])

print("Found " + str(count) + " records. Do you want to download all of them?")
download = input('[y|n]: ')
if not any(y in download for y in ["y", "Y"]):
    print('Aborting.')
    sys.exit(0)

# This fails because count is the true count of records
# while gi_list is limited by argument retmax
# assert count == len(gi_list)

# However, you also get given two additional pieces of information, the WebEnv session cookie, and the QueryKey:

webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

# Having stored these values in variables session_cookie and query_key we can use them as parameters to Bio.Entrez.efetch() instead of giving the GI numbers as identifiers.

# While for small searches you might be OK downloading everything at once, it is better to download in batches. You use the retstart and retmax parameters to specify which range of search results you want returned (starting entry using zero-based counting, and maximum number of results to return). Sometimes you will get intermittent errors from Entrez, HTTPError 5XX, we use a try except pause retry block to address this. For example,

# This assumes you have already run a search as shown above,
# and set the variables count, webenv, query_key

try:
    from urllib.error import HTTPError  # for Python 3
except ImportError:
    from urllib2 import HTTPError  # for Python 2

# For illustrative purposes, this example downloads the FASTA records in batches of three. Unless you are downloading genomes or chromosomes, you would normally pick a larger batch size.
batch_size = 100
out_handle = open(outfile_path, "w")
for start in range(0, count, batch_size):
    end = min(count, start+batch_size)
    print("Downloading record %i to %i..." % (start+1, end))
    attempt = 0
    while attempt < 3:
        attempt += 1
        try:
            fetch_handle = Entrez.efetch(db=ncbi_db,
                                         rettype=outfile_format, retmode="text",
                                         retstart=start, retmax=batch_size,
                                         webenv=webenv, query_key=query_key)
        except HTTPError as err:
            if 500 <= err.code <= 599:
                print("Received error from server %s" % err)
                print("Attempt %i of 3" % attempt)
                time.sleep(15)
            else:
                raise
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()


with open(metadata_path, mode = 'a') as f:
    f.write("User:\t" + Entrez.email + "\n" + 
            "Query Date:\t" + str(start_time) + "\n" +
            'Filename:\t' + outfile_path + "\n" +
            "Query Terms:\t" + my_query + "\n" +
            "N sequences:\t" + str(count) + "\n")
f.close()

print('Sequences written to file: ' + outfile_path)
print('Metadata written to file: ' + metadata_path)

sys.exit(0)
