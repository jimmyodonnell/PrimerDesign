# Search Field Descriptions for NCBI Sequence Databases

The details here are copied verbatim from [here](https://www.ncbi.nlm.nih.gov/books/NBK49540/), one of the NCBI's multiple, conflicting, outdated references on the subject. I have rearranged their order as they are relevant to my most common use-case, downloading sequence on the basis of the organism and gene they came from. Format is as follows:

- [Full name]  [Short name]
  - definition
  - example

### Often useful:
- [All Fields]	[ALL]
  - All terms from all search fields in the database.

- [Organism]	[ORGN]
  - The scientific and common names for the complete taxonomy of organisms that are the source of the sequence records. This vocabulary includes all available nodes in the NCBI taxonomy database.

- [Properties]	[PROP]
  - Molecular type, source database, and other properties of the sequence record. Terms indexed for this field are a useful classification system for sequence records. 
  - Molecule type
    - biomol_crna[PROP]
    - biomol_genomic[PROP]
    - biomol_mrna[PROP]
  - Cellular location
    - gene_in_genomic[PROP]
    - gene_in_mitochondrion[PROP]

- [Sequence Length]	[SLEN]
  - The total length of the sequence − the number of nucleotides or amino acids in the sequence. The colon ( : ) separates the beginning and end of a length range.
  - 100:1000[SLEN]

### Occasionally useful
- [Accession]	[ACCN]
  - The accession number assigned by NCBI.
  - AF123456[ACCN]

- [Filter]	[FILT]  [SB]
  - Filtered subsets of the database. An important kind of filter is based on the presence of links to other records. Other filters create useful subsets of data such as those set as Filters in the Discovery column of search results

- [Genome Project]
  -	The numeric unique identifier for the genome project that produced the sequence records.
  -	13139[Genome Project] (Oryza sativa Japonica)
  - 21117[Genome Project] (Pelagic Microbial Assemblages in the Oligotrophic Ocean)

- [Primary Accession]	[PACC]
  - The primary accession number of the sequence record. This is the first one appearing on the ACCESSION line in the GenBank/GenPept format. Many records have additional secondary accessions representing records that have been merged. The Accession field indexes both primary and secondary accessions.

- [Primary Organism]	[PORGN]
  - The primary organism when there is more than one source organism.

- [Publication Date]	[PDAT]
  - The date that records were made public in Entrez. The date format is YYYY/MM/DD. The colon ( : ) separates the beginning and end of a date range.
  - 2010/01:2010/12/31[PDAT]

- [SeqID String]	[SQID]
  - The NCBI identifier string for the sequence record. This is a brief structured format used by NCBI software.

- [Title]	[TI] OR [TITL]
  - Words and phrases found in the title of the sequence record. The title is the DEFINITION line of the GenBank/GenPept format of the record. This line summarizes the biology of the sequence and includes the organism, product name, gene symbol, molecule type, and sequence completeness.

### Rarely useful
- [EC/RN Number]	[ECNO]
  - Enzyme Commission (EC) number for an enzyme activity.

- [Feature Key]  [FKEY]
  - Biological features listed in the Feature Table of the sequence records. 

- [Gene Name]	[GENE]
  - Gene names annotated on database records. For NCBI Reference Sequences, these names correspond to official nomenclature guidelines when possible. Submitters provide the gene names on GenBank/GenPept records. Gene names on submitted records may be historical names or vary from official guidelines for other reasons.
  - BRCA1[GENE]

- [Keyword]	[KYWD]
  - Keywords applied by submitter or from controlled vocabularies applied by NCBI or other databases. Except for specific kinds of records, such as the examples given below, the terms in this index are not well controlled. This field is unpopulated for many GenBank/GenPept records.

- [Modification Date]	[MDAT]
  - The date of most recent modification of a sequence record. The date format is YYYY/MM/DD. Only the year is required. The Modification Date is often used as a range of dates. The colon ( : ) separates the beginning and end of a date range.

- [Molecular Weight]  [MOLWT]
  - The molecular weight in Daltons of the protein chain calculated from the amino acids only. This may not correspond to the molecular weight of the protein obtained from biological samples because of incomplete data or post-translational modifications of the protein in living systems. The colon ( : ) separates the beginning and end of a molecular weight range.

- [Protein Name]	[PROT]
  - The names of protein products as annotated on sequence records. The content of this field is not well controlled for GenBank/GenPept records and may contain inaccurate or incomplete information.

- [Substance Name]	[SUBS]
  - The names of chemical substances associated with a record. This field is only populated for sequences extracted from structure records - PDB derived sequences. The associated residue position is often included.

- [Text Word]	[WORD]
  - Text on a sequence record that is not indexed in other fields. Terms indexed here are included in an All Fields search, not generally useful.



### Publication-related
- [Author]	[AU]  [AUTH]
  - All authors from all references in the records. The format is last name [space] first initial(s), without punctuation. 
  - venter jc[AUTH]

- [Issue]	[ISS]	
  - The issue number of the journals cited on sequence records, not generally useful in sequence databases.

- [Journal]	[JOUR]
  - The name of the journals cited on sequence records. Journal names are indexed in the database in abbreviated form although many full titles are mapped to their abbreviations. Journals are also indexed by their by International Standard Serial Number (ISSN).

- [Page Number]	[PAGE]
  - The page numbers of the articles that are cited on the sequence record, not generally useful in sequence databases.

- [Volume]	[VOL]
  - Contains the volume number of the journals in references on the sequence record, not generally useful in the sequence databases.
