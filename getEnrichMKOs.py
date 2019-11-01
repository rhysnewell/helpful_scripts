#!/usr/bin/env python3
import string, sys, re, math, os, array

def main():
    try:
        genome_directory = sys.argv[1]
        output_directory = sys.argv[2]
        minimum_genes = int(sys.argv[3])
    except IndexError:
    	print("Usage: <EnrichM Gene Fastas Directory> <Output Directory> <Minimum Genes> \n",
    	"Takes output genes from EnrichM with cluster tags and rewrites each cluster to individual files \n")
    	return
    #reference_genome = sys.argv[4]
    clusters = {}
    orthologs = {}
    ko = {}
    genome_list = os.listdir(genome_directory)
    for genome_path in genome_list:
        genome = readFastaFile(genome_directory+'/'+genome_path)
        for gene in genome:
            try:
                cluster_id = gene.info.split(' ')[1]
            except IndexError:
                cluster_id = None
            try:
                ortholog_id = gene.info.split(' ')[2]
            except IndexError:
                ortholog_id = None
            try:
            	ko_id = gene.info.split(' ')[3]
            except IndexError:
            	ko_id = None
           # name = gene.name.split('_')
           # name.pop(0)
           # name.remove('contig')
           # gene.name = '_'.join(name)
            if cluster_id is not None:
                if cluster_id in clusters.keys():
            	    clusters[cluster_id].append(gene)
                else:
                    clusters[cluster_id] = [gene]
            if ortholog_id is not None:
                if ortholog_id in orthologs.keys():
                    orthologs[ortholog_id].append(gene)
                else:
                    orthologs[ortholog_id] = [gene]
            if ko_id is not None:
                if ko_id in ko.keys():
                    ko[ko_id].append([cluster_id, ortholog_id])
                else:
                    ko[ko_id] = [cluster_id, ortholog_id]
                
    try:
        os.mkdir(output_directory)
    except FileExistsError:
        print("Output directory already exists", file=sys.stderr)
        
    for cluster, genes in clusters.items():
    	if len(genes) >= minimum_genes:
            writeFastaFile(output_directory+'/'+cluster+'.fa', genes)
    for ortholog, genes in orthologs.items():
        if len(genes) >= minimum_genes:
        	writeFastaFile(output_directory+'/'+ortholog+'.fa', genes)
    with open('ko_ids.tsv', 'w') as ids:
        ids.write("cluster\tortholog\tko\n")
        for ko_id, co_ids, in ko.items():
            if type(co_ids[0]) is list:
                for co in co_ids:
                    ids.write("%s\t%s\t%s\n" % (co[0], co[1], ko_id))
            else:
                ids.write("%s\t%s\t%s\n" % (co_ids[0], co_ids[1], ko_id))
    ids.close()       

class Sequence(object):
    """ A biological sequence. Stores the sequence itself (as a compact array), 
    the alphabet (i.e., type of sequence it is), and optionally a name and further 
    information. """
    
    sequence = None # The array of symbols that make up the sequence 
    alphabet = None # The alphabet from which symbols come
    name =     None # The name (identifier) of a sequence
    info =     None # Other information (free text; e.g. annotations)
    length =   None # The number of symbols that the sequence is composed of
    gappy =    None # True if the sequence has "gaps", i.e. positions that represent deletions relative another sequence
    
    def __init__(self, sequence, alphabet = None, name = '', info = '', gappy = False):
        """ Create a sequence with the sequence data. Specifying the alphabet,
        name and other information about the sequence are all optional.
        The sequence data is immutable (stored as a string).
        Example:
        >>> myseq = Sequence('MVSAKKVPAIAMSFGVSF')
        will create a sequence with no name, and assign one of the predefined
        alphabets on the basis of what symbols were used.
        >>> myseq.alphabet.symbols
        will output the standard protein alphabet:
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
        'R', 'S', 'T', 'V', 'W', 'Y'] """
        
        self.sequence = sequence
        
        # Assign an alphabet
        # If no alphabet is provided, attempts to identify the alphabet from sequence
        self.alphabet = None
        if not alphabet is None:
            for sym in self.sequence:
                if not sym in alphabet and (sym != '-' or not gappy):  # error check: bail out
                    raise RuntimeError('Invalid symbol: %c in sequence %s' % (sym, name))
            self.alphabet = alphabet
        else:
            for alphaName in preferredOrder:
                alpha = predefAlphabets[alphaName]
                valid = True
                for sym in self.sequence:
                    if not sym in alpha and (sym != '-' or not gappy):  
                        valid = False
                        break
                if valid:
                    self.alphabet = alpha
                    break
            if self.alphabet is None:
                raise RuntimeError('Could not identify alphabet from sequence: %s' % name)
        
        # Store other information
        self.name = name
        self.info = info
        self.length = len(self.sequence)
        self.gappy = gappy
        
    def __len__(self):
        """ Defines what the "len" operator returns for an instance of Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print (len(seq))
        9
        """
        return len(self.sequence)

    def __str__(self):
        """ Defines what should be printed when the print statement is used on a Sequence instance """
        str = self.name + ': '
        for sym in self:
            str += sym
        return str
    
    def __iter__(self):
        """ Defines how a Sequence should be "iterated", i.e. what its elements are, e.g.
        >>> seq = Sequence('AGGAT', DNA_Alphabet)
        >>> for sym in seq:
                print (sym)
        will print A, G, G, A, T (each on a separate row)
        """ 
        tsyms = tuple(self.sequence)
        return tsyms.__iter__()
    
    def __contains__(self, item):
        """ Defines what is returned when the "in" operator is used on a Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print ('T' in seq)
        True
            which is equivalent to 
        >>> print (seq.__contains__('T'))
        True
        >>> print ('X' in seq)
        False
        """ 
        for sym in self.sequence:
            if sym == item:
                return True
        return False
        
    def __getitem__(self, ndx):
        """ Retrieve a specified index (or a "slice" of indices) of the sequence data.
            Calling self.__getitem__(3) is equivalent to self[3] 
        """
        if type(ndx) is slice:
            return ''.join(self.sequence[ndx])
        else:
            return self.sequence[ndx]
        
    def writeFasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        if parseDefline(self.info)[0] == self.name: # this sequence was previously "parsed" and info should hold the original header
            fasta = '>' + self.info + '\n'
        else:
            fasta = '>' + self.name + ' ' + self.info + '\n'
        data = ''.join(self.sequence)
        #nlines = int(math.ceil((len(self.sequence) - 1) / 60 + 1))
        #for i in range(nlines):
            #lineofseq = ''.join(data[i*60 : (i+1)*60]) + '\n'
        fasta += data + '\n'
        return fasta
    
    def count(self, findme = None):
        """ Get the number of occurrences of specified symbol findme OR
            if findme = None, return a dictionary of counts of all symbols in alphabet """
        if findme != None:
            cnt = 0
            for sym in self.sequence:
                if findme == sym:
                    cnt = cnt + 1
            return cnt
        else:
            symbolCounts = {}
            for symbol in self.alphabet:
                symbolCounts[symbol] = self.count(symbol)
            return symbolCounts

    def getDegapped(self):
        """ Create the sequence excluding gaps, and provide the corresponding indices for the gapped version, e.g.
        >>> gappy = Sequence('AC--TA-GA', DNA_Alphabet, name = 'myseq', gappy = True)
        >>> degapped, indices = gappy.getDegapped()
        >>> print(degapped)
            myseq: ACTAGA
        >>> print(indices)
            [0, 1, 4, 5, 7, 8]
        """
        idxs = []
        newseq = []
        for i in range(len(self.sequence)):
            if not self.sequence[i] == '-':
                newseq.append(self.sequence[i])
                idxs.append(i)
        return Sequence(newseq, self.alphabet, self.name, self.info, gappy = False), idxs

    def find(self, findme, gappy = False):
        """ Find the position of the specified symbol or sub-sequence """
        if gappy == False or self.gappy == False:
            return ''.join(self.sequence).find(findme)
        else: # if the sequence is gappy AND the function is called with gappy = True THEN run the find on the de-gapped sequence
            degapped, idxs = self.getDegapped()
            idx = ''.join(degapped).find(findme)
            return idxs[idx] if idx >= 0 else -1

"""
Below are some useful methods for loading data from strings and files.
Recognize the FASTA format (nothing fancy).
"""
def readFasta(string, alphabet = None, ignore = False, gappy = False, parse_defline = True):
    """ Read the given string as FASTA formatted data and return the list of
        sequences contained within it.
        If alphabet is specified, use it, if None (default) then guess it.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences.
        If parse_defline is False, the name will be set to everything before the first space, else parsing will be attempted."""
    seqlist = []    # list of sequences contained in the string
    seqname = None  # name of *current* sequence
    seqinfo = None
    seqdata = []    # sequence data for *current* sequence
    for line in string.splitlines():    # read every line
        if len(line) == 0:              # ignore empty lines
            continue
        if line[0] == '>':  # start of new sequence
            if seqname:     # check if we've got one current
                try:
                    current = Sequence(seqdata, alphabet, seqname, seqinfo, gappy)
                    seqlist.append(current)
                except RuntimeError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            # now collect data about the new sequence
            seqinfo = line[1:].split() # skip first char
            if len(seqinfo) > 0:
                try:
                    if parse_defline:
                        parsed = parseDefline(seqinfo[0])
                        seqname = parsed[0]
                    else:
                        seqname = seqinfo[0]
                    seqinfo = line[1:]
                except IndexError as errmsg:
                    if not ignore:
                        raise RuntimeError(errmsg)
            else:
                seqname = ''
                seqinfo = ''
            seqdata = []
        else:               # we assume this is (more) data for current
            cleanline = line.split()
            for thisline in cleanline:
                seqdata.extend(tuple(thisline.strip('*')))
    # we're done reading the file, but the last sequence remains
    if seqname:
        try:
            lastseq = Sequence(seqdata, alphabet, seqname, seqinfo, gappy)
            seqlist.append(lastseq)
        except RuntimeError as errmsg:
            if not ignore:
                raise RuntimeError(errmsg)
    return seqlist

def parseDefline(string):
    """ Parse the FASTA defline (see http://en.wikipedia.org/wiki/FASTA_format)
        GenBank, EMBL, etc                gi|gi-number|gb|accession|locus
        SWISS-PROT, TrEMBL                sp|accession|name
        ...
        Return a tuple with
        [0] primary search key, e.g. UniProt accession, Genbank GI
        [1] secondary search key, e.g. UniProt name, Genbank accession
        [2] source, e.g. 'sp' (SwissProt/UniProt), 'tr' (TrEMBL), 'gb' (Genbank)
    """
    if len(string) == 0: return ('', '', '', '')
    s = string.split()[0]
    if re.match("^sp\|[A-Z][A-Z0-9]{5}\|\S+", s):            arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^tr\|[A-Z][A-Z0-9]*\|\S+", s): arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^gi\|[0-9]*\|\S+\|\S+", s):               arg = s.split('|');  return (arg[1], arg[3], arg[0], arg[2])
    elif re.match("gb\|\S+\|\S+", s):                        arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("emb\|\S+\|\S+", s):                       arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    elif re.match("^refseq\|\S+\|\S+", s):                   arg = s.split('|');  return (arg[1], arg[2], arg[0], '')
    else: return (s, '', '', '')

def readFastaFile(filename, alphabet = None, ignore = False, gappy = False, parse_defline = True):
    """ Read the given FASTA formatted file and return the list of sequences
        contained within it. Note that if alphabet is NOT specified, it will take a
        separate guess for each sequence.
        If ignore is False, errors cause the method to fail.
        If ignore is True, errors will disregard sequence.
        If gappy is False (default), sequence cannot contain gaps,
        if True gaps are accepted and included in the resulting sequences.
        If parse_defline is False, the name will be set to everything before the first space, else parsing will be attempted."""
    fh = open(filename)
    seqlist = []
    batch = '' # a batch of rows including one or more complete FASTA entries
    rowcnt = 0
    for row in fh:
        row = row.strip()
        if len(row) > 0:
            if row.startswith('>') and rowcnt > 0:
                more = readFasta(batch, alphabet, ignore, gappy, parse_defline)
                if len(more) > 0:
                    seqlist.extend(more)
                batch = ''
                rowcnt = 0
            batch += row + '\n'
            rowcnt += 1
    if len(batch) > 0:
        more = readFasta(batch, alphabet, ignore, gappy, parse_defline)
        if len(more) > 0:
            seqlist.extend(more)
    fh.close()
    return seqlist

def writeFastaFile(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'w')
    for seq in seqs:
        fh.write(seq.writeFasta())
    fh.close()

class Alphabet(object):
    """ Defines an immutable biological alphabet (e.g. the alphabet for DNA is AGCT)
    that can be used to create sequences (see sequence.py).
    We use alphabets to define "tuple" tables, where entries are keyed by combinations
    of symbols of an alphabet (see class TupleStore below).
    Alphabets are used to define probability distributions for stochastic events
    (see prob.py). """

    def __init__(self, symbolString):
        """ Construct an alphabet from a string of symbols. Lower case characters
        will be converted to upper case, repeated characters are ignored.
        Example of constructing the DNA alphabet:
        >>> alpha = Alphabet('ACGTttga')
        >>> alpha.symbols
        ('A', 'C', 'G', 'T') """

        # Add each symbol to the symbols list, one at a time, and ignore doubles (could use "set" here...)
        _symbols = [] # create a temporary list
        for s in symbolString:
            if not str(s).upper()[0] in _symbols:
                _symbols.append(str(s).upper()[0])
        _symbols.sort() # we put them in alphabetical (one canonical) order
        # OK done extracting, put them in place
        self.symbols = tuple(_symbols); # create the immutable tuple from the extracted list
        self.length = len(self.symbols)
        self.annotations = {}

    def __str__(self):
        return str(self.symbols)

    def __len__(self):
        return len(self.symbols)

    def __iter__(self):
        return self.symbols.__iter__()

    def __getitem__(self, ndx):
        """ Retrieve the symbol(s) at the specified index (or slice of indices) """
        return self.symbols[ndx]

    def __contains__(self, sym):
        """ Check if the given symbol is a member of the alphabet. """
        return sym in self.symbols

    def index(self, sym):
        """ Retrieve the index of the given symbol in the alphabet. """
        # If the symbol is valid, use the tuple's index function
        if sym in self.symbols:
            syms = self.symbols
            return syms.index(sym)
        else:
            raise RuntimeError('Symbol %s is not indexed by alphabet %s' % (sym, str(self.symbols)))

    def __eq__(self, rhs):
        """ Test if the rhs alphabet is equal to ours. """
        if rhs == None:
            return False
        if len(rhs) != len(self):
            return False
        # OK we know they're same size...
        for sym in self.symbols:
            if not sym in rhs:
                return False
        return True

    def isSubsetOf(self, alpha2):
        """ Test if this alphabet is a subset of alpha2. """
        for sym in self.symbols:
            if not sym in alpha2:
                return False
        return True

    def isSupersetOf(self, alpha2):
        """ Test if this alphabet is a superset of alpha2. """
        return alpha2.isSubsetOf(self)

    def annotateSym(self, label, sym, value):
        try:
            lookup = self.annotations[label]
        except KeyError:
            lookup = self.annotations[label] = {}
        lookup[sym] = value

    def annotateAll(self, label, symdictOrFilename):
        if isinstance(symdictOrFilename, str): # we assume it is a filename
            fh = open(symdictOrFilename)
            string = fh.read()
            d = {}
            for line in string.splitlines():
                if len(line.strip()) == 0:
                    continue
                sections = line.split()
                symstr, value = sections[0:2]
                for sym in symstr:
                    d[sym] = value
            fh.close()
        else: # we assume it is a dictionary
            d = symdictOrFilename
        for sym in d:
            self.annotateSym(label, sym, d[sym])

    def getAnnotation(self, label, sym):
        try:
            lookup = self.annotations[label]
            return lookup[sym]
        except KeyError:
            return None


""" Below we declare alphabets that are going to be available when
this module is imported """
Bool_Alphabet = Alphabet('TF')
DNA_Alphabet = Alphabet('ACGT')
DNA_Alphabet_wN = Alphabet('ACGTN')
RNA_Alphabet_wN = Alphabet('ACGUN')
RNA_Alphabet = Alphabet('ACGU')
Protein_Alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
Protein_Alphabet_wX = Protein_wX = Alphabet('ACDEFGHIKLMNPQRSTVWYX')
Protein_Alphabet_wSTOP = Protein_wSTOP = Alphabet('ACDEFGHIKLMNPQRSTVWY*')
Protein_wGAP  = Alphabet('ACDEFGHIKLMNPQRSTVWY-')
DSSP_Alphabet = Alphabet('GHITEBSC')
DSSP3_Alphabet = Alphabet('HEC')

predefAlphabets = {'Bool_Alphabet': Bool_Alphabet,
                   'DNA': DNA_Alphabet,
                   'RNA': RNA_Alphabet,
                   'DNAwN': RNA_Alphabet_wN,
                   'RNAwN': DNA_Alphabet_wN,
                   'Protein': Protein_Alphabet,
                   'ProteinwX': Protein_wX,
                   'ProteinwSTOP' : Protein_wSTOP,
                   'ProteinwGAP': Protein_wGAP,
                   'DSSP_Alphabet' : DSSP_Alphabet,
                   'DSSP3_Alphabet' : DSSP3_Alphabet}
# The preferred order in which a predefined alphabet is assigned to a sequence
# (e.g., we'd want to assign DNA to 'AGCT', even though Protein is also valid)
preferredOrder = ['Bool_Alphabet', 'DNA', 'RNA', 'DNAwN', 'RNAwN', 'Protein', 'ProteinwX', 'ProteinwSTOP',
                  'ProteinwGAP', 'DSSP_Alphabet', 'DSSP3_Alphabet']
# Useful annotations
DNA_Alphabet.annotateAll('html-color', {'A':'green','C':'orange','G':'red','T':'#66bbff'})
RNA_Alphabet.annotateAll('html-color', {'A':'green','C':'orange','G':'red','U':'#66bbff'})

Protein_Alphabet.annotateAll('html-color', {'G':'green', 'C': 'green', 'P':'green','S':'#66bbff','T': '#66bbff','H':'red',
                                            'K':'red','R':'red','F':'green','Y':'#66bbff','W':'green','I':'green',
                                            'L':'green','M':'green', 'N': '#66bbff', 'Q': '#66bbff', 'V':'green',
                                            'A' : 'green', 'D': 'orange', 'E': 'orange'})


# ------------------ Substitution Matrix ------------------

class TupleStore(dict):
    """ Internal utility class that can be used for associating
    a value with ordered n-tuples (n=1..N).
    Read/write functions are defined for instances of this class.
    """

    def __init__(self, alphas=None, entries=None, sparse=True):
        """
        Manage entries keyed by symbol-tuples with values of arbitrary type.
        If alphas is None, the alphabet(s) are inferred from the provided entries.
        If entries is None, all entries are defined by possible combinations of symbols from specified alphabets,
        and are assumed to be None until specified. Either alphas or entries must be supplied.
        If sparse is True, a sparse memory-saving encoding is used, if false, a time-saving, more flexible encoding is used.
        >>> matrix = TupleStore({'AA': 2, 'AW': -3, 'WW': 4, 'AR': -1})
        >>> matrix[('A', 'W')]
        -3
        >>> matrix['AR']
        -1
        """
        assert sparse, "Currently only sparse encoding is implemented."
        assert alphas or entries, "Either alphabets or entries (from which alphabets can be inferred) must be supplied."
        self.sparse = sparse         # sparse encoding if true
        if alphas == None:
            self.alphas = None       # need to figure out alphabet from supplied entries
            self.keylen = None       # tuple length not known yet
        elif type(alphas) is Alphabet:
            self.alphas = tuple ([ alphas ]) # make it into a tuple
            self.keylen = 1          # tuple length 1
        else:
            self.alphas = alphas     # alphabets are supplied
            self.keylen = len(alphas)# length of tuples is the same as the number alphabets

        # Check if entries are supplied to the constructor
        if entries == None:
            self.entries = entries = {}
        elif type(entries) is dict:
            raise RuntimeError("When specified, entries must be a dictionary")
        # Check length of tuples, must be the same for all
        for entry in entries:
            if self.keylen == None:
                self.keylen = len(entry)
            elif self.keylen != len(entry):
                raise RuntimeError("All entries must have the same number of symbols")

        # go through each position in tuples, to check what alphabet is right
        myalphas = []                   # my suggestions from entries (need to be subsets of specified)
        for idx in range(self.keylen):
            symset = set()              # we collect all symbols in position idx here
            for key in entries:
                symset.add(key[idx])
            myalpha = Alphabet(symset)
            myalphas.append(myalpha)
            if self.alphas != None:     # if specified it needs to be a superset of that we constructed
                if not self.alphas[idx].isSupersetOf(myalpha):
                    raise RuntimeError("Specified alphabet is not compatible with specified entries")

        if self.alphas == None:     # if not specified to constructor use those we found
            self.alphas = tuple(myalphas)

        for key in entries:
            self[key] = entries[key]

    def _isValid(self, symkey):
        for idx in range(self.keylen):
            if not symkey[idx] in self.alphas[idx]:
                return False
        return True

    def __setitem__(self, symkey, value):
        assert self.keylen == len(symkey), "All entries in dictionary must be equally long"
        assert self._isValid(symkey), "Invalid symbol in entry"
        self.entries[symkey] = value

    def __getitem__(self, symkey):
        """ Return the score matching the given symbols together."""
        assert self.keylen == len(symkey), "Entries must be of the same length"
        try:
            return self.entries[symkey]
        except KeyError:
            return None

    def __iadd__(self, symkey, ivalue):
        assert self.keylen == len(symkey), "All entries in dictionary must be equally long"
        assert self._isValid(symkey), "Invalid symbol in entry"
        try:
            self.entries[symkey] += ivalue
        except KeyError:
            self.entries[symkey] = ivalue

    def __isub__(self, symkey, ivalue):
        assert self.keylen == len(symkey), "All entries in dictionary must be equally long"
        assert self._isValid(symkey), "Invalid symbol in entry"
        try:
            self.entries[symkey] -= ivalue
        except KeyError:
            self.entries[symkey] = -ivalue

    def getAll(self, symkey=None):
        """ Return the values matching the given symbols together.
        symkey: tuple (or list) of symbols or None (symcount symbol); if tuple is None, all entries are iterated over.
        """
        if symkey == None:
            symkey = []
            for idx in range(self.keylen):
                symkey.append(None)
        else:
            assert self.keylen == len(symkey), "Entries must be of the same length"
        for idx in range(self.keylen):
            if symkey[idx] != None:
                if not symkey[idx] in self.alphas[idx]:
                    raise RuntimeError("Invalid entry: must be symbols from specified alphabet or None")
        return TupleEntries(self, symkey)

    def __iter__(self):
        return TupleEntries(self, tuple([None for _ in range(self.keylen)]))

    def items(self, sort = False):
        """ In a dictionary-like way return all entries as a list of 2-tuples (key, prob).
        If sort is True, entries are sorted in descending order of value.
        Note that this function should NOT be used for big (>5 variables) tables."""
        ret = []
        for s in self.entries:
            if self[s] != None:
                ret.append((s, self[s]))
        if sort:
            return sorted(ret, key=lambda v: v[1], reverse=True)
        return ret

class TupleEntries(object):
    """ Iterator class for multiple entries in a tuple store.
    """
    def __init__(self, tuplestore, symkey):
        self.tuplestore = tuplestore
        self.symkey = symkey
        self.symcount = []
        self.indices = []
        for ndx in range(tuplestore.keylen):
            if symkey[ndx] == None:
                self.indices.append(ndx)
                self.symcount.append(0)        # start at this index to alter symbol
            else:
                self.symcount.append(None)     # do not alter this symbol
        self.nextIsLast = False

    def __iter__(self):
        return self

    def __next__(self):
        """ Step through sequence of entries, either
        (if not sparse) with a step-size based on alphabet-sizes and what symbols are specified or
        (if sparse) with calls to tuple store based on all possible symbol combinations."""

        if self.nextIsLast:
            raise StopIteration

        mykey = [] # construct current combination from known and unspecified symbols
        for ndx in range(self.tuplestore.keylen):
            if (self.symkey[ndx] == None):
                sym = self.tuplestore.alphas[ndx][self.symcount[ndx]]
                mykey.append(sym)
            else:
                mykey.append(self.symkey[ndx])

        # decide which ndx that should be increased (only one)
        self.nextIsLast = True # assume this is the last round (all counters are re-set)
        for ndx in self.indices:
            if self.symcount[ndx] == len(self.tuplestore.alphas[ndx]) - 1: # if we just entered the last symbol of this alphabet
                self.symcount[ndx] = 0                  # reset count here
            else:
                self.symcount[ndx] = self.symcount[ndx] + 1
                self.nextIsLast = False
                break

        return tuple(mykey)

if __name__ == '__main__':
	main()