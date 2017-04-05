#!/usr/bin/env python

def getUniqueSeqs(inReads, outFasta):
    uniqueSeqs = {}
    uniqueReads = []
    outFile = open(outFasta, 'w')
    for read in inReads:
        if read.seq not in uniqueSeqs:
            uniqueSeqs[read.seq] = [read]
        else:
            uniqueSeqs[read.seq].append(read)
    for u in sorted(uniqueSeqs, key=lambda k: len(uniqueSeqs[k]), reverse=True):
        firstRead = uniqueSeqs[u][0]
        outFile.write(firstRead.header + ';size=' + str(len(uniqueSeqs[u])) + \
                        ';\n' + firstRead.seq + '\n')
    outFile.close()
    return uniqueSeqs

#Translate unoise2 output and SINTAX output into headers of original reads
#Input: reads imported from denoised otus
#	taxonomic dictionary imported from sintax output
#	unique dictionary created by getUniqueSeqs
#Output: a file is written with primer-filtered reads, plus otu and tax
#	information in their headers
def otuToHeaders(otuReads, taxInfo, unique, outFileName):
    outFile = open(outFileName, 'w')
    for otu in otuReads:
        for read in unique[otu.seq]:
            otuID = otu.header.replace('>','').split(';')[0]
            finalTax = taxInfo[otu.header.replace('>','')][1]
            newHeader = read.header + ';' + otuID + ';tax=' + finalTax + ';'
            outFile.write(newHeader + '\n')
            outFile.write(read.seq + '\n')
    outFile.close()

#Written specifically for uparse output
def uniqueSeqsToOTU(uparseMap):
    otuToHeaders = {}
    i = 1
    with open(uparseMap, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if line[1] == 'OTU':
                otuToHeaders['OTU' + str(i)] = [line[0]]
                i += 1
    with open(uparseMap, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if (line[1] == 'match') and ('top=OTU' in line[2]):
                otuID = line[2].split('top=')[1].split('(')[0].strip()
                otuToHeaders[otuID].append(line[0])
    return otuToHeaders

#Translate unoise2 output and SINTAX output into headers of original reads
#Input: reads imported from denoised otus
#	taxonomic dictionary imported from sintax output
#	unique dictionary created by getUniqueSeqs
#Output: a file is written with primer-filtered reads, plus otu and tax
#	information in their headers
def otuToHeaders(otuReads, taxInfo, unique, outFileName):
    outFile = open(outFileName, 'w')
    for otu in otuReads:
        for read in unique[otu.seq]:
            otuID = otu.header.replace('>','').split(';')[0]
            finalTax = taxInfo[otu.header.replace('>','')][1]
            newHeader = read.header + ';' + otuID + ';tax=' + finalTax + ';'
            outFile.write(newHeader + '\n')
            outFile.write(read.seq + '\n')
    outFile.close()
