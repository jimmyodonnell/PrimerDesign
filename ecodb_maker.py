#!/usr/bin/env python

'''
--------------------------------------------------------------------------------
Original program written by Riaz et al.?

Attempted improvements by:
Author: Jimmy O'Donnell <jodonnellbio@gmail.com>

NOTES: 
  - removed option to use taxdb; I could find no documentation about its usage
  - removed progress bar; it didn't work right. Use module progressbar2

TODO: if no taxonomy, download and unzip
TODO: write metadata file:
  - input filename:
  - date:
  - for each sequence, write organisms and length

'''

import re
import gzip
import struct
import sys
import argparse
import os
from datetime import datetime

parser = argparse.ArgumentParser(
  description = 'Create a database for use with ecoPrimers and ecoPCR')

parser.add_argument('-i', '--input', 
  help = 'input file containing DNA sequences in GenBank format', 
  required = True)

parser.add_argument('-t', '--taxonomy',
  help = 'path to folder containing NCBI taxonomy files', 
  required = True)

parser.add_argument('-o', '--output',
  help = 'path to output', 
  required = True)

def infile_type_checker(infilepath):
    
    with open(infilepath, 'r') as f:
        first_line = f.readline()
        first_part = str.split(first_line)[0]

    # guess_method = "contents"
    if first_part == 'LOCUS':
        file_type = 'genbank'
    elif first_part == 'ID':
        file_type = 'embl'
    elif first_part[0] == '>':
        file_type = 'fasta'
    else:
        # guess_method = "extension"
        infile_ext = os.path.splitext(infilepath)[1].lower()
        if infile_ext in ('.gb'):
            file_type = 'genbank'
        if infile_ext in ('.embl'):
            file_type = 'embl'
        if infile_ext in ('.fa', '.fsa', '.fasta'):
            file_type = 'fasta'
        else:
            raise ValueError('could not guess input file type')
    
    return(file_type)

#####
#
#
# Generic file function
#
#
#####

def universalOpen(file):
    if isinstance(file,str):
        if file[-3:] == '.gz':
            rep = gzip.open(file)
        else:
            rep = open(file)
    else:
        rep = file
    return rep

def universalTell(file):
    if isinstance(file, gzip.GzipFile):
        file=file.myfileobj
    return file.tell()


#####
#
#
# NCBI Dump Taxonomy reader
#
#
#####

def endLessIterator(endedlist):
    for x in endedlist:
        yield x
    while(1):
        yield endedlist[-1]
   
class ColumnFile(object):
    
    def __init__(self,stream,sep=None,strip=True,types=None):
        if isinstance(stream,str):
            self._stream = open(stream)
        elif hasattr(stream,'next'):
            self._stream = stream
        else:
            raise ValueError,'stream must be string or an iterator'
        self._delimiter=sep
        self._strip=strip
        if types:
            self._types=[x for x in types]
            for i in xrange(len(self._types)):
                if self._types[i] is bool:
                    self._types[i]=ColumnFile.str2bool
        else:
            self._types=None
            
    def str2bool(x):
        return bool(eval(x.strip()[0].upper(),{'T':True,'V':True,'F':False}))
                    
    str2bool = staticmethod(str2bool)
            
        
    def __iter__(self):
        return self
    
    def next(self):
        line = self._stream.next()
        data = line.split(self._delimiter)
        if self._strip or self._types:
            data = [x.strip() for x in data]
        if self._types:
            it = endLessIterator(self._types)
            data = [x[1](x[0]) for x in ((y,it.next()) for y in data)]
        return data
    
def taxonCmp(t1,t2):
    if t1[0] < t2[0]:
        return -1
    elif t1[0] > t2[0]:
        return +1
    return 0

def bsearchTaxon(taxonomy,taxid):
    taxCount = len(taxonomy)
    begin = 0
    end   = taxCount 
    oldcheck=taxCount
    check = begin + end / 2
    while check != oldcheck and taxonomy[check][0]!=taxid :
        if taxonomy[check][0] < taxid:
            begin=check
        else:
            end=check
        oldcheck=check
        check = (begin + end) / 2
        
        
    if taxonomy[check][0]==taxid:
        return check
    else:
        return None
        
    
    
def readNodeTable(file):

    file = universalOpen(file)
    
    nodes = ColumnFile(file, 
                       sep='|', 
                       types=(int,int,str,
                              str,str,bool,
                              int,bool,int,
                              bool,bool,bool,str))
    sys.stderr.write("Reading taxonomy dump file...\n")
    taxonomy=[[n[0],n[2],n[1]] for n in nodes]
    sys.stderr.write("List all taxonomy rank...\n")    
    ranks =list(set(x[1] for x in taxonomy))
    ranks.sort()
    ranks = dict(map(None,ranks,xrange(len(ranks))))
    
    sys.stderr.write("Sorting taxons...\n")
    taxonomy.sort(taxonCmp)

    sys.stderr.write("Indexing taxonomy...\n")
    index = {}
    for t in taxonomy:
        index[t[0]]=bsearchTaxon(taxonomy, t[0])
    
    sys.stderr.write("Indexing parent and rank...\n")
    for t in taxonomy:
        t[1]=ranks[t[1]]
        t[2]=index[t[2]]
        
        
    return taxonomy,ranks,index

def nameIterator(file):
    file = universalOpen(file)
    names = ColumnFile(file, 
                       sep='|', 
                       types=(int,str,
                              str,str))
    for taxid,name,unique,classname,white in names:
        yield taxid,name,classname
        
def mergedNodeIterator(file):
    file = universalOpen(file)
    merged = ColumnFile(file, 
                       sep='|', 
                       types=(int,int,str))
    for taxid,current,white in merged:
            yield taxid,current
    
def deletedNodeIterator(file):
    file = universalOpen(file)
    deleted = ColumnFile(file, 
                       sep='|', 
                       types=(int,str))
    for taxid,white in deleted:
            yield taxid
    
def readTaxonomyDump(taxdir):
    taxonomy,ranks,index = readNodeTable('%s/nodes.dmp' % taxdir)
    
    sys.stderr.write("Adding scientific name...\n")

    alternativeName=[]
    for taxid,name,classname in nameIterator('%s/names.dmp' % taxdir):
        alternativeName.append((name,classname,index[taxid]))
        if classname == 'scientific name':
            taxonomy[index[taxid]].append(name)
        
    sys.stderr.write("Adding taxid alias...\n")
    for taxid,current in mergedNodeIterator('%s/merged.dmp' % taxdir):
        index[taxid]=index[current]
    
    sys.stderr.write("Adding deleted taxid...\n")
    for taxid in deletedNodeIterator('%s/delnodes.dmp' % taxdir):
        index[taxid]=None
    
    return taxonomy,ranks,alternativeName,index

#####
#
#
#  Genbank/EMBL sequence reader
#
#
#####

def entryIterator(file):
    file = universalOpen(file)
    rep =[]
    for line in file:
        rep.append(line)
        if line == '//\n':
            rep = ''.join(rep)
            yield rep
            rep = []
            
def fastaEntryIterator(file):
    file = universalOpen(file)
    rep =[]
    for line in file:
        if line[0] == '>' and rep:
            rep = ''.join(rep)
            yield rep
            rep = []
        rep.append(line)
    if rep:
        rep = ''.join(rep)
        yield rep
    
_cleanSeq = re.compile('[ \n0-9]+')
            
def cleanSeq(seq):
    return _cleanSeq.sub('',seq)
    
    
_gbParseID = re.compile('(?<=^LOCUS {7})[^ ]+(?= )',re.MULTILINE)   
_gbParseDE = re.compile('(?<=^DEFINITION {2}).+?\. *$(?=[^ ])',re.MULTILINE+re.DOTALL)   
_gbParseSQ = re.compile('(?<=^ORIGIN).+?(?=^//$)',re.MULTILINE+re.DOTALL)  
_gbParseTX = re.compile('(?<= /db_xref="taxon:)[0-9]+(?=")')
  
def genbankEntryParser(entry):
    Id = _gbParseID.findall(entry)[0]
    De = ' '.join(_gbParseDE.findall(entry)[0].split())
    Sq = cleanSeq(_gbParseSQ.findall(entry)[0].upper())
    try:
        Tx = int(_gbParseTX.findall(entry)[0])
    except IndexError:
        Tx = None
    return {'id':Id,'taxid':Tx,'definition':De,'sequence':Sq}

######################

_cleanDef = re.compile('[\nDE]')

def cleanDef(definition):
    return _cleanDef.sub('',definition)

_emblParseID = re.compile('(?<=^ID {3})[^ ]+(?=;)',re.MULTILINE)   
_emblParseDE = re.compile('(?<=^DE {3}).+?\. *$(?=[^ ])',re.MULTILINE+re.DOTALL)   
_emblParseSQ = re.compile('(?<=^  ).+?(?=^//$)',re.MULTILINE+re.DOTALL)  
_emblParseTX = re.compile('(?<= /db_xref="taxon:)[0-9]+(?=")')

def emblEntryParser(entry):
    Id = _emblParseID.findall(entry)[0]
    De = ' '.join(cleanDef(_emblParseDE.findall(entry)[0]).split())
    Sq = cleanSeq(_emblParseSQ.findall(entry)[0].upper())
    try:
        Tx = int(_emblParseTX.findall(entry)[0])
    except IndexError:
        Tx = None
    return {'id':Id,'taxid':Tx,'definition':De,'sequence':Sq}


######################

_fastaSplit=re.compile(';\W*')

def parseFasta(seq):
    seq=seq.split('\n')
    title = seq[0].strip()[1:].split(None,1)
    id=title[0]
    if len(title) == 2:
        field = _fastaSplit.split(title[1])
    else:
        field=[]
    info = dict(x.split('=',1) for x in field if '=' in x)
    definition = ' '.join([x for x in field if '=' not in x])
    seq=(''.join([x.strip() for x in seq[1:]])).upper()   
    return id,seq,definition,info

  
def fastaEntryParser(entry):
    id,seq,definition,info = parseFasta(entry)
    Tx = info.get('taxid',None)   
    if Tx is not None:
        Tx=int(Tx)
    return {'id':id,'taxid':Tx,'definition':definition,'sequence':seq}
    

def sequenceIteratorFactory(entryParser,entryIterator):
    def sequenceIterator(file):
        for entry in entryIterator(file):
            yield entryParser(entry)
    return sequenceIterator


def taxonomyInfo(entry,connection):
    taxid = entry['taxid']
    curseur = connection.cursor()
    curseur.execute("""
                        select taxid,species,genus,family,
                               taxonomy.scientificName(taxid) as sn,
                               taxonomy.scientificName(species) as species_sn,
                               taxonomy.scientificName(genus) as genus_sn,
                               taxonomy.scientificName(family) as family_sn
                        from
                            (   
                             select alias                      as taxid,
                               taxonomy.getSpecies(alias) as species,
                               taxonomy.getGenus(alias)   as genus,
                               taxonomy.getFamily(alias)  as family
                                from taxonomy.aliases
                               where id=%d ) as tax
                    """ % taxid)
    rep = curseur.fetchone()
    entry['current_taxid']=rep[0]
    entry['species']=rep[1]
    entry['genus']=rep[2]
    entry['family']=rep[3]
    entry['scientific_name']=rep[4]
    entry['species_sn']=rep[5]
    entry['genus_sn']=rep[6]
    entry['family_sn']=rep[7]
    return entry
    
#####
#
#
# Binary writer
#
#
#####
    
def ecoSeqPacker(sq):
    
    compactseq = gzip.zlib.compress(sq['sequence'],9)
    cptseqlength  = len(compactseq)
    delength   = len(sq['definition'])
    
    totalSize = 4 + 20 + 4 + 4 + 4 + cptseqlength + delength
    
    packed = struct.pack('> I I 20s I I I %ds %ds' % (delength,cptseqlength),
                         totalSize,
                         sq['taxid'],
                         sq['id'],
                         delength,
                         len(sq['sequence']),
                         cptseqlength,
                         sq['definition'],
                         compactseq)
    
    assert len(packed) == totalSize+4, "error in sequence packing"
    
    return packed

def ecoTaxPacker(tx):
    
    namelength = len(tx[3])
    
    totalSize = 4 + 4 + 4 + 4 + namelength
    
    packed = struct.pack('> I I I I I %ds' % namelength, 
                         totalSize, 
                         tx[0],
                         tx[1],
                         tx[2], 
                         namelength,
                         tx[3])
    
    return packed

def ecoRankPacker(rank):
    
    namelength = len(rank)
    
    packed = struct.pack('> I %ds' % namelength,
                         namelength,
                         rank)
    
    return packed
                
def ecoNamePacker(name):
    
    namelength = len(name[0])
    classlength= len(name[1])
    totalSize =  namelength + classlength + 4 + 4 + 4 + 4
    
    packed = struct.pack('> I I I I I %ds %ds' % (namelength,classlength),
                         totalSize,
                         int(name[1]=='scientific name'),
                         namelength,
                         classlength,
                         name[2],
                         name[0],
                         name[1])
    
    return packed
    
def ecoSeqWriter(file,input,taxindex,parser):
    output = open(file,'wb')
    input  = universalOpen(input)
    entries = parser(input)
    seqcount=0
    skipped = []

    output.write(struct.pack('> I',seqcount))
    
    for entry in entries:
        if entry['taxid'] is not None:
            try:
                entry['taxid']=taxindex[entry['taxid']]
            except KeyError:
                entry['taxid']=None
            if entry['taxid'] is not None:
                seqcount+=1
                output.write(ecoSeqPacker(entry))
            else:
                skipped.append(entry['id'])
            where = universalTell(input)
            # progressBar(where, inputsize)
            # sys.stderr.write(" Readed sequences : %d     " % seqcount,
        else:
            skipped.append(entry['id'])
        
    output.seek(0,0)
    output.write(struct.pack('> I',seqcount))
    
    output.close()
    return skipped
        

def ecoTaxWriter(file,taxonomy):
    output = open(file,'wb')
    output.write(struct.pack('> I',len(taxonomy)))
    
    for tx in taxonomy:
        output.write(ecoTaxPacker(tx))

    output.close()
    
def ecoRankWriter(file,ranks):
    output = open(file,'wb')
    output.write(struct.pack('> I',len(ranks)))

    rankNames = ranks.keys()
    rankNames.sort()
    
    for rank in rankNames:
        output.write(ecoRankPacker(rank))

    output.close()

def nameCmp(n1,n2):
    name1=n1[0].upper()
    name2=n2[0].upper()
    if name1 < name2:
        return -1
    elif name1 > name2:
        return 1
    return 0


def ecoNameWriter(file,names):
    output = open(file,'wb')
    output.write(struct.pack('> I',len(names)))

    names.sort(nameCmp)
    
    for name in names:
        output.write(ecoNamePacker(name))

    output.close()

def infile_parser_picker(infile_type):
    if infile_type == 'genbank':
        infile_parser = sequenceIteratorFactory(genbankEntryParser, 
                            entryIterator)
    elif infile_type == 'fasta':
        infile_parser = sequenceIteratorFactory(fastaEntryParser, 
                            fastaEntryIterator)
    elif infile_type == 'embl':
        infile_parser = sequenceIteratorFactory(emblEntryParser, 
                            entryIterator)
    else:
        raise ValueError("could not pick an infile parser")
    return(infile_parser)

def ecoDBWriter(prefix,taxonomy,seqFileNames):
    
    sys.stderr.write("Writing database...\n")
    
    ecoRankWriter('%s.rdx' % prefix, taxonomy[1])
    ecoTaxWriter('%s.tdx' % prefix, taxonomy[0])
    ecoNameWriter('%s.ndx' % prefix, taxonomy[2])
  
    filecount = 0
    for filename in seqFileNames:
        infile_type = infile_type_checker(filename)
        the_parser = infile_parser_picker(infile_type)
        filecount+=1
        sk=ecoSeqWriter('%s_%03d.sdx' % (prefix,filecount), 
                     filename, 
                     taxonomy[3], 
                     the_parser)
        if sk:
            sys.stderr.write("Skipped entry :\n")
            sys.stderr.write(sk + "\n")
        
if __name__ == '__main__':
    
    full_time_start = datetime.now()

    args = vars(parser.parse_args())

    taxon_time_start = datetime.now()

    taxonomy = readTaxonomyDump(args['taxonomy'])
    
    taxon_time_end = datetime.now()

    duration_taxonomy = str(taxon_time_end - taxon_time_start)

    sys.stderr.write("Taxonomy compilation completed in: " + duration_taxonomy + "\n")
    
    infile_type = infile_type_checker(args['input'])
    the_parser = infile_parser_picker(infile_type)
    ecoDBWriter(args['output'], taxonomy, [args['input']])

    full_time_end = datetime.now()
    
    duration_full = str(full_time_end - full_time_start)
    
    sys.stderr.write("Completed in: " + duration_full + "\n")
