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
    currentSeq = ''
    currentHeader = ''
    currentCluster = ''
    for line in inputFile:
	if '>' in line:
	    currentSeq = line.split(' ')[0].replace('>','')
	    currentHeader = line.strip()
	    currentCluster = line.split(' ')[1].split('#')[0]
	else:
	    readObj = fastaSeq(currentSeq, currentCluster, currentHeader, line.strip())
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
	newRead = fastaSeq(read.seq_name, read.cluster, read.header,
			splitSeqRev[0])
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
	newRead = fastaSeq(read.seq_name, read.cluster, read.header,
			   splitSeq[1])
	usableReads.append(newRead)
    return usableReads

# Trim sequences to a specified length
def trimLength(reads, outputLength):
    returnReads = []
    for read in reads:
	if len(read.seq) >= outputLength:
	    newSeq = read.seq[0:outputLength]
	    newRead = fastaSeq(read.seq_name, read.cluster, read.header, newSeq)
	    returnReads.append(newRead)
    return returnReads

# Select reads based on provided sample list, return selected reads
def selectReads(sampList, reads):
    outReads = []
    for s in sampList:
        for r in reads:
            if s in r.seq_name:
                outReads.append(r)
    return outReads

class fastaSeq(object):
    def __init__(self, seq_name, cluster, header, seq):
	self.seq_name = seq_name
	self.cluster = cluster
	self.header = header
	self.seq = seq
