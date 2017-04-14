# Primer Design

You will need:

1. R, including the package taxize 
2. python, including the module Bio ('BioPython')
3. [ecoPrimers](https://git.metabarcoding.org/obitools/ecoprimers/wikis/home); Documentation [here](http://metabarcoding.org/obitools/doc/scripts/ecoPrimers.html)
4. [ecoPCR](https://git.metabarcoding.org/obitools/ecopcr/wikis/home); documentation [here](http://metabarcoding.org/obitools/doc/scripts/ecoPCR.html)
5. NCBI Taxonomy database dump file [available here](ftp://ftp.ncbi.nih.gov/pub/taxonomy). Or, from the command line:

`wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz`

Your life will be easier if you get `ecoPrimers`, `ecoPCR`, and `ecoPCRFormat.py` into your shell's `$PATH` variable. Accomplishing that is [up to you](https://www.google.com/search?q=how+to+add+to+path+variable).

## Here are the steps:

1. Download sequences in GenBank format that include target and non-target taxon. (python)
2. Convert those into a format that ecoPrimers and ecoPCR can recognize (~~`ecoPCRFormat.py`~~ `ecodb_maker.py`)
3. Identify primers that will work on the target and not the others (`ecoPrimers`, `taxdump`)
4. Refine primer list (R).
5. Evaluate primer performance (ecoPCR)
6. ???
7. Profit

### 1. Download sequences from GenBank

#### `get_genbank.py`

You will need a collection of DNA sequences *in GenBank format* (*.gb) from both the target and non-target taxa. You could access the GenBank interface using a web browser, and accomplish this task. The downside of this approach is that you could easily forget the details of the sequences you gathered: what the search query was, when it was downloaded, et cetera. If you ~~are a glutton for punishment~~ want to keep track of these details, use a script. 

However, downloading sequences using a script or the command line is not a straightforward task, for several reasons. First, NCBI views the process as two very discrete steps: querying their databases and retrieving data. Second, the system NCBI requires you to use ('Entrez') is not terribly straightforward both for constructing your search query (i.e. the stuff you type into a search box) and for accessing the data.
Third, NCBI puts the onus on users to not overload their servers, and will restrict your access if you violate their terms. This means you have to do a little fiddling to avoid being banned.
Finally, the documentation for the whole process, like all NCBI documentation, is shit.

The least painful way to go about this is in python, using the module `Bio` (aka "BioPython").
I recall BioPython installation being confusing when I first did it a few years back. 
Now, with a little more experience, it was relatively painless -- I installed the latest version of it for Python 3 using `pip3 install numpy; pip3 install biopython`.

Once you've got biopython installed, the script `get_genbank.py` should be relatively easy to use. It requires two arguments: your email address, and your search query.
The queries look weird (e.g. `'"Homo sapiens"[Organism] AND BRCA1[GENE] NOT 21117[Genome Project]'`), and I haven't put time into making the query construction any easier.
You can see all of the query fields in the file `ncbi_query_fields.md`, where you might notice there are many possible fields to search.
These include some that are redundant ([WORD] and [ALL]), and others which could be misleading (e.g. [PROT]: *"The content of this field is not well controlled for GenBank/GenPept records and may contain inaccurate or incomplete information."*).

Remember that NCBI's sequence database (a.k.a. GenBank a.k.a. `nt`) is not meticulously curated.
If your query was `'"Oophaga pumilio"[ORGN] AND MT-CO1[GENE]'`, you might expect to get all of the sequences from the mitochondrial gene encoding cytochrome c oxidase subunit 1 in the taxon *Oophaga pumilio*.
However, you'd miss any instances of that sequence where it is part of a whole mitochondrial genome, or where the submitter simply used an incorrect gene ID, or none at all.
Thus, you'd do better to use `'"Oophaga pumilio"[ORGN] AND ("cytochrome c oxidase subunit"[ALL] OR CO1[ALL] OR COI[ALL] OR COX1[ALL] OR MT-CO1[ALL])'` or `'"Oophaga pumilio"[ORGN] AND cytochrome[ALL] NOT "cytochrome b"[ALL]'`.

`¯\_(ツ)_/¯`

### 2. Convert sequences into another format

#### ~~`ecoPCRFormat.py`~~ `ecodb_maker.py`

**tl;dr:** I modified `ecoPCRFormat.py` to make it easier to use. I call my version `ecodb_maker.py`. Usage is simple, and it prints help if you need it. Use it like so:  

`./ecodb_maker.py -i sequences.gb -t /path/to/taxdump -o output_name`)

If you thought that was more painful than it should be, [hold on to your butts](https://www.youtube.com/watch?v=-W6as8oVcuM). 
A team in France that does a lot of great metabarcoding work has written a set of programs that are widely used, but not terribly user friendly. 
Their usage, inner workings, and output are poorly documented.
This could be improved by incorporating these programs as, for example, a python module. 
They might already be part of the python module `obitools`, but again, poor documentation leaves me not really understanding how they are related to one another or to this package.

`ecoPrimers` requires that your sequences be formatted in a specific way. 
I would have expected this to happen as part of the program; however, they require the user to first run a separate script to create the database: `ecoFormat.py`. 
I suppose the reasoning is that formatting a database could take a while for lots of sequences, or having multiple redundant databases is an inefficient use of disk space? 
In any case, it is what it is, and you have to run it.
According to the help for `ecoPCRFormat.py`, it should accept FASTA format input, but I could only get it to accept input in GenBank format. 

Formatting a database in any sort of logical manner is beyond me. 
It's a pretty standard thing in software to give a computer program directions on where it can find input and to where it should direct its output.
Yet, as best I can tell, this is not the logic under which `ecoPCRFormat.py` operates. 
Here are some of the problems I ran into.

In the following examples, I'm using shell variables for the input file (`$GB_FILE`), the `taxdump` folder (`$TAX`), and the relative path to the output files (`$OUT`).
I tried a variety of output files. 
The path to the output file should be a basename with no extension.
If you do not specify an output name, the output files will have a basename the same as the input file (including extension) but with an additional extension appended.
If you specify a directory that does not exist, it will fail.

The argument `--taxonomy` must point to the folder called 'taxdump', not any of the files within it.
Rather than check all of the arguments first, `ecoPCRFormat.py` starts by reading the entirety of the NCBI taxonomy dump.
When it appears something is successfully happening with the taxonomy database, the script outputs the following lines (I replace them with `*** builds taxonomy ***` in my examples).

```
Reading taxonomy dump file...
List all taxonomy rank...
Sorting taxons...
Indexing taxonomy...
Indexing parent and rank...
Adding scientific name...
Adding taxid alias...
Adding deleted taxid...
```

Here are some of my misadventures:

Simply running the script with no options gives an obscure error message, rather than providing help or example usage: 
```
$ ecoPCRFormat.py
Traceback (most recent call last):
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 644, in <module>
    if opt['taxmod']=='dump':
KeyError: 'taxmod'
```

The same error message is reported for the following screw ups:

```
$ ecoPCRFormat.py --genbank $OUT --taxonomy $TAX --name  $GB_FILE
$ ecoPCRFormat.py --genbank $GB_FILE --name $OUT --taxonomy $TAX
```

The argument for output is "--name" (and can't be a directory, especially one that doesn't exist), and the input file is given simply as the last argument. Thus, the following attempts all bomb:

```
$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX $OUT
*** builds taxonomy ***
Traceback (most recent call last):
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 650, in <module>
    ecoDBWriter(opt['prefix'], taxonomy, filenames, opt['parser'])
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 574, in ecoDBWriter
    parser)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 486, in ecoSeqWriter
    input  = universalOpen(input)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 29, in universalOpen
    rep = open(file)
IOError: [Errno 2] No such file or directory: 'lolwut'
```

I then made a directory using the name of $OUT

```
$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX $OUT
*** builds taxonomy ***
Traceback (most recent call last):
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 650, in <module>
    ecoDBWriter(opt['prefix'], taxonomy, filenames, opt['parser'])
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 574, in ecoDBWriter
    parser)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 486, in ecoSeqWriter
    input  = universalOpen(input)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 29, in universalOpen
    rep = open(file)
IOError: [Errno 21] Is a directory: 'lolwut'
```

And then tried pointing to a subdirectory within it where I wanted the output...

```
$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX lolwut/lolwut
*** builds taxonomy ***
Traceback (most recent call last):
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 650, in <module>
    ecoDBWriter(opt['prefix'], taxonomy, filenames, opt['parser'])
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 574, in ecoDBWriter
    parser)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 486, in ecoSeqWriter
    input  = universalOpen(input)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 29, in universalOpen
    rep = open(file)
IOError: [Errno 2] No such file or directory: 'lolwut/lolwut'

$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX lolwut*
*** builds taxonomy ***
Traceback (most recent call last):
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 650, in <module>
    ecoDBWriter(opt['prefix'], taxonomy, filenames, opt['parser'])
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 574, in ecoDBWriter
    parser)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 486, in ecoSeqWriter
    input  = universalOpen(input)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 29, in universalOpen
    rep = open(file)
IOError: [Errno 21] Is a directory: 'lolwut'

$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX lolwut/lol*
*** builds taxonomy ***
Traceback (most recent call last):
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 650, in <module>
    ecoDBWriter(opt['prefix'], taxonomy, filenames, opt['parser'])
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 574, in ecoDBWriter
    parser)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 486, in ecoSeqWriter
    input  = universalOpen(input)
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 29, in universalOpen
    rep = open(file)
IOError: [Errno 2] No such file or directory: 'lolwut/lol*'
```

I'm not sure why this failed...
```
$ ecoPCRFormat.py --genbank -t "${TAX}" --name "lolwut/lol" "${GB_FILE}" 
*** builds taxonomy ***
Traceback (most recent call last):
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 650, in <module>
    ecoDBWriter(opt['prefix'], taxonomy, filenames, opt['parser'])
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 564, in ecoDBWriter
    ecoRankWriter('%s.rdx' % prefix, taxonomy[1])
  File "/Users/threeprime/bin/ecoPCRFormat.py", line 530, in ecoRankWriter
    output = open(file,'wb')
IOError: [Errno 2] No such file or directory: 'lolwut/lol.rdx'
```

... but I finally realized the final argument should be the input file.
The following attempts do not produce an error:

```
$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX 
*** builds taxonomy ***

$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX > $OUT
*** builds taxonomy ***
```

It output a *.ndx, *.rdx, and *.tdx* file; but successful runs appear to also generate a `*.sdx` file, which this does not. So apparently, you can create a valid `ecoP*` database without any sequence data whatsoever?
No error message, though.

The following *does* work:
```
$ ecoPCRFormat.py --genbank --taxonomy $TAX --name $OUT $GB_FILE

$ ecoPCRFormat.py --genbank --name $GB_FILE --taxonomy $TAX /Users/threeprime/Desktop/gbdl_170301-141036*
```

A successful run looked something like this after building the taxonomy files:
```
1.5 % |-                            ] remain : 00:00:00 Readed sequences : 1  
1.5 % |-                            ] remain : 00:00:00 Readed sequences : 2  
3.3 % |#-                           ] remain : 00:00:00 Readed sequences : 3  
3.3 % |#-                           ] remain : 00:00:00 Readed sequences : 4  
3.3 % |#-                           ] remain : 00:00:00 Readed sequences : 5  
3.3 % |#-                           ] remain : 00:00:00 Readed sequences : 6  
5.1 % |##-                          ] remain : 00:00:00 Readed sequences : 7  
5.1 % |##-                          ] remain : 00:00:00 Readed sequences : 8  
6.9 % |###-                         ] remain : 00:00:00 Readed sequences : 9  
6.9 % |###-                         ] remain : 00:00:00 Readed sequences : 1  
6.9 % |###-                         ] remain : 00:00:00 Readed sequences : 1  
8.7 % |####-                        ] remain : 00:00:00 Readed sequences : 1  
8.7 % |####-                        ] remain : 00:00:00 Readed sequences : 1  
...
94.2 % |#########################-   ] remain : 00:00:00 Readed sequences : 1 
99.6 % |###########################- ] remain : 00:00:00 Readed sequences : 1
100.0 % |############################/] remain : 00:00:00 Readed sequences : 114     

```


Here are some things I learned:
1. Your input file is given as the last argument without a flag.
2. The input options (`--genbank`, `--embl`, `--fasta`) simply set the type -- do not give the input file as a value to these arguments.
3. If no output name is given, the output will be the name of your input (including the extension) plus the new extensions.
4. You can give the output a new name, including in a subdirectory, but the directory must exist.
5. The output is much larger than the input. When input was 563 KB, output was 212400 KB (~377x).
6. Output files are binary, not text, making it difficult to summarize them in a straightforward way.


### 3. Find primers that might work

#### `ecoPrimers`

Once you have formatted a database, you can run `ecoPrimers`. 
This program generates sequences that meet some given criteria given your input sequences -- e.g. they are of a given length range, they match target taxa but exclude non-targets, they are a given distance apart, etc. 
For some reason, it relies on you having a local copy of NCBI's taxonomy database dump, `taxdump`. 
I imagine this is to identify which clades are nested within others.
Similarly, because the structure of this database identifies taxa as *numbers* rather than *names*, you'll need to have handy the NCBI Taxon ID for both the target taxon and the taxon you'd like to exclude. 
The simple script `ncbi_taxid_get.r` can do this for you, if you provide it with a scientific name.

Several of the options are confusing.
The following options tripped me up for longer than I'd like to admit: 
```
-r <TAXID>
Defines the example sequence records (example dataset). Only the sequences of the corresponding taxonomic group identified by its TAXID are taken into account for designing the barcodes and the primers. The TAXID is an integer that can be found either in the NCBI taxonomic database, or using the ecofind program.

-i <TAXID>
Defines the counterexample sequence records (counterexample dataset). The barcodes and primers will be selected in order to avoid the counterexample taxonomic group identified by its TAXID.

-E <TAXID>
Defines an counterexample taxonomic group (identified by its TAXID) within the example dataset.
```
Can you spot the subtle difference between the second (`-i`) and third (`-E`) options? 
I couldn't. 
`-r` specifies the ingroup or target taxon.
I think the source of confusion is that they refer to sequences from the target/ingroup taxon as the 'example dataset', whereas I think of the example dataset as the entire database I assembled to run the program.
An example of the intended use case might be that you generate a database (using `ecodb_maker.py`) containing every mitochondrial genome from genbank, and then you want to design primers that only amplify vertebrates (`-r`), exclude arthropods (`-i`), and also exclude rodents (`-E`).
Keep in mind that the whole database could also contain birds, fish, etc. 
Of course, I could be way off.


Additionally, there is an argument to turn on "double strand mode" (`-D`), and a separate argument to turn on "single strand mode" (`-S`). 
Presumably, then, one could run it in both double strand and single strand mode? Seems odd.

The output of `ecoPrimers` is a tab-separated table providing the sequences and some data about each one. 
This table can be pretty large, perhaps because my input parameters were too lenient. 
Thus, before you go buying primers and testing them in the lab or even in silico, you should further refine the criteria for good primers.


### 3. Select best primers

#### `ecoprimer_filter.r`

Not much to explain here; just use some basic filters to whittle the output of `ecoPrimers` down to only the primers that are most likely to work.
This is where you select for characteristics of of the primers not covered by ecoPrimers that you might , including:

- Primers are between 18 and 30 bp long
- Amplicon between 70 and 150 base pairs
- GC content 40-60%
- Repeats of no more than 4 base pairs
- Melting temperature between 50 and 65 ºC
- Melting temperature no different than 2 ºC
- Ensure no 3' complementarity that may generate primer-dimer

- Picking a probe: 
  - Melting temp ~8-10 ºC higher than primers
  - shorter than 30 bp
  - no 5' G

There are lots of useful, but conflicting sources of information on these considerations ([here is one](https://www.idtdna.com/pages/decoded/decoded-articles/pipet-tips/decoded/2011/09/12/steps-for-a-successful-qpcr-experiment)).

### 4. Test the primers in silico

#### `ecoPCR`

I have no idea what this does differently than `ecoPrimers`.


### 5. ???

After running `ecoPrimers`, you should engage in the art of black magic known as 'lab work'.


### 6. Profit

Self-explanatory. Don't forget to step on all the little people along the way.

