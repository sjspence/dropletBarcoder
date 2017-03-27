from optparse import OptionParser
import re

degenerate = {'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT',
		'M':'AC', 'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG',
		'N':'ACGT'}
alt_map = {'ins':'0'}
complement = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'Y':'R', 'R':'Y',
		'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'V':'B', 
		'D':'H', 'H':'D', 'N':'N'} 

def revComplement(seq):    
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def makeRE(sequence):
    exportSeq = ''
    for base in sequence:
        if base not in degenerate:
            exportSeq += base
        else:
            d = degenerate[base]
            d = '[' + d + ']'
            exportSeq += d
    return exportSeq

def getHeaderParams(header):
    spaceSplit = header.split(' ')
    currentSeq = spaceSplit[0].strip().replace('>','')
    currentHeader = header.strip()
    if (len(spaceSplit) > 1) and ('#' in spaceSplit[1]):
        currentCluster = spaceSplit[1].split('#')[0]
    else:
        currentCluster = None
    return currentHeader, currentSeq, currentCluster

# Import data from fasta file and save in header and sequence variables
def importFasta(inputFileName):
    reads = []
    inputFile = open(inputFileName,'r')
    currentDNA = ''
    currentSeq = ''
    currentHeader = ''
    currentCluster = ''
    for line in inputFile:
        if '>' in line:
            if currentDNA != '':
                readObj = FastaSeq(currentHeader, currentDNA, currentSeq,
					currentCluster)
                reads.append(readObj)
            currentHeader, currentSeq, currentCluster = getHeaderParams(line)
            currentDNA = ''
        else:
            currentDNA += line.strip()
    if currentDNA != '':
        readObj = FastaSeq(currentHeader, currentDNA, currentSeq,
				currentCluster)
        reads.append(readObj)
    inputFile.close()
    return reads

# Parse the joined fasta reads for designed primer sequence structure
# Output: sequences that match the designed structure
def removeFwdRevPrimer(reads, fwd, rev):
    usableReads = []
    forward = makeRE(fwd)
    reverse = makeRE(rev)
    for read in reads:
        if re.search(forward, read.seq[0:len(forward)]):
            pass
        else:
            continue
        splitSeq = re.split(forward, read.seq)
        fwdRemoved = splitSeq[1]
        if re.search(reverse,
                read.seq[len(read.seq)-len(reverse):len(read.seq)]):
            pass
        else:
            continue
        splitSeqRev = re.split(reverse,fwdRemoved)
        newRead = FastaSeq(read.header, splitSeqRev[0], read.seq_id, 
                                read.cluster)
        usableReads.append(newRead)
    return usableReads

# Parse the joined fasta reads for designed primer sequence structure
# Output: sequences that match the designed structure
def removeFwdPrimer(reads, fwd):
    usableReads = []
    primer = makeRE(fwd)
    for read in reads:
        if not re.search(primer, read.seq[0:len(primer)]):
            continue
        splitSeq = re.split(primer, read.seq)
        if len(splitSeq[1]) == 0:
            continue
        newRead = FastaSeq(read.header, splitSeq[1], read.seq_id,
                                read.cluster)
        usableReads.append(newRead)
    return usableReads

# Parse the joined fasta reads for designed primer sequence structure
# Output: sequences that match the designed structure
def filtBarcodePrimers(reads, bcLength, fwd, rev):
    usableReads = []
    forward = makeRE(fwd)
    reverse = makeRE(revComplement(rev))
    endFwdLoc = bcLength + len(forward)
    for read in reads:
        if (re.search(forward, read.seq[bcLength:endFwdLoc]) and 
		re.search(reverse,
		read.seq[len(read.seq)-len(reverse):len(read.seq)])):
            pass
        else:
            continue
        fwdRemoved = re.split(forward, read.seq)[1]
        splitSeqRev = re.split(reverse, fwdRemoved)[0]
        newHeader = read.header + ' droplet_bc=' + read.seq[0:bcLength]
        newRead = FastaSeq(newHeader, splitSeqRev, read.seq_id,
                                read.cluster)
        usableReads.append(newRead)
    return usableReads

# Split a read by a provided degenerate primer sequence
def splitByDegenerate(seq, degenerateSeq):
    REseq = makeRE(degenerateSeq)
    newSeq = re.split(REseq, seq)
    return newSeq

# Trim sequences to a specified length
def trimLength(reads, outputLength):
    returnReads = []
    for read in reads:
        if len(read.seq) >= outputLength:
            newSeq = read.seq[0:outputLength]
            newRead = FastaSeq(read.header, newSeq, read.seq_id,
                                read.cluster)
            returnReads.append(newRead)
    return returnReads

# Select reads based on provided sample list, return selected reads
def selectSamples(sampList, reads):
    outReads = []
    for s in sampList:
        for r in reads:
            if s == r.seq_id.split('_')[0]:
                outReads.append(r)
    return outReads

# Remove reads that map to a provided sample list, return remaining reads
def removeSamples(sampList, reads):
    outReads = []
    for r in reads:
        match = False
        for s in sampList:
            if s == r.seq_id.split('_')[0]:
                match = True
        if not match:
            outReads.append(r)
    return outReads

class FastaSeq(object):
    #def __init__(self, seq_name=None, cluster=None, header, seq):
    def __init__(self, header, seq, seq_id=None, cluster=None):
        self.header = header
        self.seq = seq
        self.seq_id = seq_id
        self.cluster = cluster
