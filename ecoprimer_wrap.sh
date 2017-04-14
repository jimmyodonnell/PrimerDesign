#!/usr/bin/env bash

# usage: ./ecoprimer_wrap.sh "/path/to/ecoprimer_db_files_no_extension" "/path/to/output_no_extension"

# also run: ./src/ecoPrimers -h

# source: http://metabarcoding.org/obitools/doc/scripts/ecoPrimers.html

# OPTIONS

#===============================================================================
# REQUIRED:
#-------------------------------------------------------------------------------
# Path to ecoP- formatted database, EXCLUDING EXTENSION OR NUMBER
ECO_DB="${1}"
# -d <filename>
# Filename containing the reference sequence records used for designing the barcode markers and primers (see obiconvert for a description of the database format).
#===============================================================================

# command Ryan used:
# ecoprimers -d ecodb -l 100 -L 500 -e 2 -r 50794 -E 100191 -t species -s 0.95 -x 0.1 -q 0.95 -O 22 -c -3 5

#===============================================================================
# OPTIONAL
ERROR_MAX="2" # RPK:2
# -e <INTEGER>
# Maximum number of errors (mismatches) allowed per primer (default: 0).

LEN_INSERT_MIN="30" # RPK:100
# -l <INTEGER>
# Minimum length of the barcode, excluding primers.

LEN_INSERT_MAX="150" # RPK:500
# -L <INTEGER>
# Maximum length of the barcode, excluding primers.

INGROUP_TAXID="50794" # RPK:Cymatogaster aggregata(50794); Embiotocide:50791
# -r <TAXID>
# Defines the example sequence records (example dataset). Only the sequences of the corresponding taxonomic group identified by its TAXID are taken into account for designing the barcodes and the primers. The TAXID is an integer that can be found either in the NCBI taxonomic database, or using the ecofind program.

OUTGROUP_TAXID="0" # 0 to ignore 
# -i <TAXID>
# Defines the counterexample sequence records (counterexample dataset). The barcodes and primers will be selected in order to avoid the counterexample taxonomic group identified by its TAXID.

INGROUP_EXCLUDE_TAXID="100191" # RPK:Hyperprosopon(100191) # 0 for no excluded in-taxon
# -E <TAXID>
# Defines an counterexample taxonomic group (identified by its TAXID) within the example dataset.

CIRCULAR="YES" # RPK: YES
# -c
# Considers that the sequences of the database are circular (e.g. mitochondrial or chloroplast DNA).

STRICT_MATCH_3P="5" # RPK:5
# -3 <INTEGER>
# Defines the number of nucleotides on the 3’ end of the primers that must have a strict match with their target sequences.

INGROUP_MATCH_PROP="0.95" # RPK:0.95
# -q <FLOAT>
# Defines the strict matching quorum, i.e. the proportion of the sequence records in which a strict match between the primers and their targets occurs (default: 0.7)

INGROUP_PARAM_PROP="0.95" # RPK:0.95
# -s <FLOAT>
# Defines the sensitivity quorum, i.e. the proportion of the example sequence records that must fulfill the specified parameters for designing the barcodes and the primers.

OUTGROUP_MATCH_PROP="0.1" # RPK:0.1
# -x <FLOAT>
# Defines the false positive quorum, i.e. the maximum proportion of the counterexample sequence records that fulfill the specified parameters for designing the barcodes and the primers.

TAXON_LEVEL="species" # RPK:species
# -t <TAXONOMIC_LEVEL>
# Defines the taxonomic level that is considered for evaluating the barcodes and primers in the output of ecoPrimers. The default taxonomic level is the species level. When using a taxonomic database builts from a NCBI taxonomy dump files, the other possible taxonomic levels are genus, family, order, class, phylum, kingdom, and superkingdom.

DOUBLE_STRANDED="NO" # FIXME Why are there two arguments for this?
# -D
# Sets the double strand mode.

SINGLE_STRANDED="NO" # FIXME Why are there two arguments for this?
# -S
# Sets the single strand mode.

PRIMER_LENGTH="22" # RPK:22
# -O <INTEGER>
# Sets the primer length (default: 18).

MELT_TEMP_METHOD="1"
# -m <1|2>
# Defines the method used for estimating the Tm (melting temperature) between the primers and their corresponding target sequences (default: 1).
# 1 SantaLucia method (SantaLucia J (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS, 95, 1460-1465).
# 2 Owczarzy method (Owczarzy R, Vallone PM, Gallo FJ et al. (1997) Predicting sequence-dependent melting stability of short duplex DNA oligomers. Biopolymers, 44, 217-239).

SALT_CONC="0.05"
# -a <FLOAT>
# Salt concentration used for estimating the Tm (default: 0.05).

MULTIMATCH="NO"
# -U
# No multi match of a primer on the same sequence record.

REFERENCE_SEQ="0"
# -R <TEXT>
# Defines the reference sequence by indicating its identifier in the database.

#=====================================================================
# IRRELEVANT FOR SCRIPTING
# -A
# Prints the list of all identifiers of sequence records in the database.

# -f
# Remove data mining step during strict primer identification.

# -v
# Stores statistic file about memory usage during strict primer identification.

# -h
# Print help.
#===============================================================================

START_TIME=$(date +%Y%m%d_%H%M)

if [ "${#2}" = 0 ]; then
  OUTBASE="ecoprimers_""$START_TIME"
else
  OUTBASE="${2}"
fi

OUTFILE="$OUTBASE"".tsv"
METADATA="$OUTBASE"".md"

ARGS=" -d $ECO_DB \
  -e "${ERROR_MAX}" \
  -l "${LEN_INSERT_MIN}" \
  -L "${LEN_INSERT_MAX}" \
  -r "${INGROUP_TAXID}" \
  -3 "${STRICT_MATCH_3P}" \
  -q "${INGROUP_MATCH_PROP}" \
  -s "${INGROUP_PARAM_PROP}" \
  -x "${OUTGROUP_MATCH_PROP}" \
  -t "${TAXON_LEVEL}" \
  -O "${PRIMER_LENGTH}" \
  -m "${MELT_TEMP_METHOD}" \
  -a "${SALT_CONC}" 
"

# ADD FLAGS
if [[ "${OUTGROUP_TAXID}" -gt 0 ]]; then
  ARGS="${ARGS} -i ${OUTGROUP_TAXID}"
fi

if [[ "${INGROUP_EXCLUDE_TAXID}" -gt 0 ]]; then
  ARGS="${ARGS} -E ${INGROUP_EXCLUDE_TAXID}"
fi

if [[ "${CIRCULAR}" == "YES" ]]; then
  ARGS="${ARGS} -c "
fi

if [[ "${DOUBLE_STRANDED}" == "YES" ]]; then
  ARGS="${ARGS} -D "
fi

if [[ "${SINGLE_STRANDED}" == "YES" ]]; then
  ARGS="${ARGS} -S "
fi

if [[ "${MULTIMATCH}" == "YES" ]]; then
  ARGS="${ARGS} -U "
fi

if [[ "${REFERENCE_SEQ}" != "0" ]]; then
  ARGS="${ARGS} -R ${REFERENCE_SEQ} "
fi

ecoPrimers ${ARGS} > "${OUTFILE}"


# OUTPUT FORMAT
# The output file contains several columns, with ‘|’ as separator, and describes the characteristics of each barcode and its associated primers.

# column 1: serial number
# column 2: sequence of primer 1
# column 3: sequence of primer 2
# column 4: Tm (melting temperature) of primer 1, without mismatch
# column 5: lowest Tm of primer 1 against example sequence records
# column 6: Tm of primer 2, without mismatch
# column 7: lowest Tm of primer 2 against example sequence records
# column 8: number of C or G in primer 1
# column 9: number of C or G in primer 2
# column 10: GG (Good-Good) means that both primer are specific to the example dataset, GB or BG (Good-Bad or Bad-Good) means that only one of the two primers is specific to the example dataset
# column 11: number of sequence records of the example dataset that are properly amplified according to the specified parameters
# column 12: proportion of sequence records of the example dataset that are properly amplified according to the specified parameters
# column 13: yule-like output
# column 14: number of taxa of the example dataset that are properly amplified according to the specified parameters
# column 15: number of taxa of the counterexample dataset that are properly amplified according to the specified parameters
# column 16: proportion of taxa of the example dataset that are properly amplified according to the specified parameters (Bc index)
# column 17: number of taxa of the example dataset that are properly identified
# column 18: proportion of taxa of the example dataset that are properly identified (Bs index)
# column 19: minimum length of the barcode in base pairs for the example sequence records (excluding primers)
# column 20: maximum length of the barcode in base pairs for the example sequence records (excluding primers)
# column 21: average length of the barcode in base pairs for the example sequence records(excluding primers)


# ecoPCRFormat:::

# 1. make directory for output with timestamp

# 2. make the ecodb/ecodb.*dx files
# ecoPCRFormat ...

# 3. Write a metadata file that includes the sequences IDS:
# ecoPrimers -d /Users/threeprime/Desktop/lolwut/lol -A 

# YAML:
# --- # Indented Block
#  name: John Smith
#  age: 33
# input filename:
# date:

echo "SUMMARY:
Using the program ecoPrimers (CITE ECOPRIMERS), against a database of HOWMANY sequences of WHICHTAXA WHICHLOCI, we designed primers which target NCBITAXID ""${INGROUP_TAXID}"" but not ""${INGROUP_EXCLUDE_TAXID}"" or ""${OUTGROUP_TAXID}"". 
We required primers strictly matched ""${INGROUP_MATCH_PROP}"" of the target taxon's sequences, and that ""${INGROUP_PARAM_PROP}"" of the target taxon's sequences must fulfill other primer design parameters (evaluated at the level of ""${TAXON_LEVEL}""). 
Primers were excluded if they matched more than ""${OUTGROUP_MATCH_PROP}"" of the non-target sequences in the database. 

The primers were designed to amplify a fragment between ""${LEN_INSERT_MIN}"" and ""${LEN_INSERT_MAX}"" base pairs in size excluding primers (""${PRIMER_LENGTH}"" base pairs each), with a maximum of ""${ERROR_MAX}"" errors (mismatches) allowed per primer. 
We required a strict match of at least ""${STRICT_MATCH_3P}"" base pairs on the 3' ends. 
Melting temperature was evaluated using method CHECKTHISAGAINSTECOPRIMERSHELP ""${MELT_TEMP_METHOD}"" at a salt concentration of ""${SALT_CONC}"".

The full command was:
ecoprimers "$ARGS"
" > "${METADATA}"

exit 0
