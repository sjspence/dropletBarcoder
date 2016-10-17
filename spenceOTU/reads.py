from optparse import OptionParser
import re

degenerate = {'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT',
		'M':'AC', 'B':'CGT', 'D':'AGT', 'H':'ACT', 'V':'ACG',
		'N':'ACGT'}

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
	    spaceSplit = line.split(' ')
	    currentSeq = spaceSplit[0].replace('>','')
	    currentHeader = line.strip()
	    if (len(spaceSplit) > 1) and ('#' in spaceSplit[1]):
		currentCluster = spaceSplit[1].split('#')[0]
	    else:
		currentCluster = None
	    currentDNA = ''
	else:
	    if '>' in next(inputFile):
		currentDNA += line.strip()
		readObj = fastaSeq(currentHeader, currentDNA, currentSeq, currentCluster)
		reads.append(readObj)
	    else:
		currentDNA += line.strip()
    inputFile.close()
    return reads

def importFasta(inputFileName):
    reads = []
    inputFile = open(inputFileName, 'r')
    currentSeq = ''
    currentHeader = ''
    currentCluster = ''
    for line in inputFile:


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
	newRead = fastaSeq(read.header, splitSeqRev[0], read.seq_id, 
				read.cluster)
	usableReads.append(newRead)
    return usableReads

# Parse the joined fasta reads for designed primer sequence structure
# Output: sequences that match the designed structure
def removeFwdPrimer(reads, fwd):
    usableReads = []
    primer = makeRE(fwd)
    for read in reads:
        if re.search(primer, read.seq[0:len(primer)]):
            pass
        else:
            continue
        splitSeq = re.split(primer, read.seq)
	newRead = fastaSeq(read.header, splitSeq[1], read.seq_id,
				read.cluster)
	usableReads.append(newRead)
    return usableReads

# Trim sequences to a specified length
def trimLength(reads, outputLength):
    returnReads = []
    for read in reads:
	if len(read.seq) >= outputLength:
	    newSeq = read.seq[0:outputLength]
	    newRead = fastaSeq(read.header, newSeq, read.seq_id,
				read.cluster)
	    returnReads.append(newRead)
    return returnReads

# Select reads based on provided sample list, return selected reads
def selectReads(sampList, reads):
    outReads = []
    for s in sampList:
        for r in reads:
            if s in r.seq_id:
                outReads.append(r)
    return outReads

class fastaSeq(object):
    #def __init__(self, seq_name=None, cluster=None, header, seq):
    def __init__(self, header, seq, seq_id=None, cluster=None):
	self.header = header
	self.seq = seq
	self.seq_id = seq_id
	self.cluster = cluster
