alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def revComplement(seq):    
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def pairConcatenate(reads1, reads2):
    pairedReads = []
    clusters1 = {}
    clusters2 = {}
    for r in reads1:
        clusters1[r.cluster] = r
    for r in reads2:
        clusters2[r.cluster] = r
    for c in clusters1:
        if c in clusters2:
            newSeq = clusters1[c].seq + revComplement(clusters2[c].seq)
            newSeqN = clusters1[c].seq + 'N' + revComplement(clusters2[c].seq)
            newRead = pairedSeq(clusters1[c].seq_id, clusters2[c].seq_id,
                                c, clusters1[c].header, clusters2[c].header,
                                newSeq, newSeqN)
            pairedReads.append(newRead)
    return pairedReads

def buildClusterDict(reads):
    clusterDict = {}
    for r in reads:
        clusterDict[r.cluster] = r
    return clusterDict

#Takes two lists of reads and returns a list of [read1, read2] read objects from the same sequencing cluster
def readPairList(reads1, reads2):
    pairs = []
    clusters1 = buildClusterDict(reads1)
    clusters2 = buildClusterDict(reads2)
    for c in clusters1:
        if c in clusters2:
            pairs.append([clusters1[c], clusters2[c]])
    return pairs

class pairedSeq(object):
    def __init__(self, seq_name1, seq_name2, cluster, header1, header2, seq, 
		seqN):
        self.seq_name1 = seq_name1
        self.seq_name2 = seq_name2
        self.cluster = cluster
        self.header1 = header1
        self.header2 = header2
        self.seq = seq
        self.seqN = seqN
