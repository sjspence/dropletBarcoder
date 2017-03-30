#!/usr/bin/env python

import reads as rd
import taxonomy as tx

def makeTreeConstraint(inFileName, outFileName):
    inFile = open(inFileName, 'r')
    taxAssignments = []
    for line in inFile:
        x = line.strip().split('\t')
        taxAssignments.append(x)
    inFile.close()
    seqID = []
    tax = []
    for i in taxAssignments:
        seqID.append(i[0])
        tax.append(i[1].split(';'))
    taxSet = [[] for i in range(7)]    #kingdom, phylum, class, order, family, genus, species
    for a in tax:
        for i in range(7):
            if len(a) > i and a[i] not in taxSet[i]:
                taxSet[i].append(a[i])
    #now make a list of each binary identifier
    outFile = open(outFileName, 'w')
    test = []
    posNum = []
    for x in taxAssignments:
        taxBinary = [[] for i in range(7)]
        pos = taxAssignments.index(x)
        test.append(seqID[pos])
        posNum.append(pos)
        for i in range(7):
            for j in taxSet[i]:
                k = "0"
                if j in x[1]:
                    k = "1"
                taxBinary[i].append(k)
        tBins=[]
        for w in taxBinary:
            tBins.append(''.join(w))
        outFile.write(">" + seqID[pos] + "\n" + ''.join(tBins) + "\n")
    outFile.close()

def alignmentToSequence(inFileName, outFileName):
    inFile = open(inFileName, 'r')
    bases = {'A':'A', 'U':'T', 'T':'T', 'G':'G', 'C':'C', 'N':'N', '\n':'',
		'\r':'', '-':'', '.':''}
    header = ''
    sequence = ''
    for line in inFile:
        if '>' in line:
            header = line
        else:
            for i in line:
                sequence += bases[i]
    outFile = open(outFileName, 'w')
    outFile.write(header + sequence)
    inFile.close()
    outFile.close()

#Choose representative sequence based on abundance in taxonomic group, export fasta of representative seqs
#INPUT:  usearch denoised fasta file of biological sequences
#        usearch sintax file of taxonomic assignments
#        file path to write representative sequences to
#OUTPUT: representative tOTU sequences written to specified file
def tOTU_pickRepSeqs(denoisedFile, sintaxFile, outFile):
    denoised = rd.importFasta(denoisedFile)        #REMOVE EBS IN NEXT THREE LINES
    taxDict = tx.importSintax(sintaxFile, 'final')
    otuDf = tx.tOTUmap(taxDict) 
    taxSeqs = {}
    for read in denoised:
        counts = int(read.header.split('size=')[1].replace(';',''))
        zOTU = read.seq_id.split(';')[0]
        tOTU = otuDf['tOTU'][zOTU]
        if tOTU not in taxSeqs:
            taxSeqs[tOTU] = [[read, counts]]
        else:
            taxSeqs[tOTU].append([read, counts])
    outReads = []
    for t in taxSeqs:
        maxCt = 0
        maxRead = ''
        for i in taxSeqs[t]:
            if i[1] > maxCt:
                maxCt = i[1]
                maxRead = i[0]
        #Edit header to include tOTU designation
        maxRead.header = maxRead.header.replace('>', '>' + t + ' ')
        outReads.append(maxRead)
    rd.exportFasta(outReads, outFile)
