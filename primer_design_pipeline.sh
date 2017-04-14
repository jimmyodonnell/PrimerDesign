#!/usr/bin/env bash

################################################################################
# Author: Jimmy O'Donnell <jodonnellbio@gmail.com>
################################################################################

# set a basename for files
FILEBASE="embio_wa_cytb"

#-------------------------------------------------------------------------------
# Download sequences from GenBank
#-------------------------------------------------------------------------------
MYQUERY='"Embiotocidae"[organism] NOT Ditrema[organism] NOT Neoditrema[organism] AND gene_in_mitochondrion[PROP] AND ("cytochrome b") NOT ("COI")'

./get_genbank.py -e jimmyod@uw.edu -q "${MYQUERY}" -f gb -o "${FILEBASE}"

#-------------------------------------------------------------------------------
# Create eco* database
#-------------------------------------------------------------------------------
TAXDUMP="/Users/jimmy.odonnell/Data/NCBI/databases/taxonomy/taxdump"
./ecodb_maker.py --input "${FILEBASE}".gb -t "$TAXDUMP" -o "${FILEBASE}"


#-------------------------------------------------------------------------------
# Create an alignment
#-------------------------------------------------------------------------------
# make a fasta copy
./gb_to_fasta.py -i "${FILEBASE}".gb

# correct orientation while aligning
mafft --adjustdirectionaccurately --maxiterate 1000 --localpair \
  "${FILEBASE}".fasta > "${FILEBASE}".aln

#-------------------------------------------------------------------------------
# Run ecoprimers
#-------------------------------------------------------------------------------
./ecoprimer_wrap.sh "${FILEBASE}" "${FILEBASE}"_ecoprimers

#-------------------------------------------------------------------------------
# Filter ecoprimer results
#-------------------------------------------------------------------------------
Rscript ecoprimer_filter.R "${FILEBASE}"_ecoprimers.tsv

#-------------------------------------------------------------------------------
# Align primers to database
#-------------------------------------------------------------------------------
mafft --addfragments "${FILEBASE}"_ecoprimers_primers.fasta --thread -1 \
  --adjustdirectionaccurately "${FILEBASE}".aln > "${FILEBASE}"_primers.aln

#-------------------------------------------------------------------------------
# Run ecopcr
#-------------------------------------------------------------------------------
ecopcr -d "${FILEBASE}" -l 30 -L 150 -e 2 -k -c

