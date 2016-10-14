def exportFasta(reads, outFileName):
    outFile = open(outFileName, 'w')
    for read in reads:
	outFile.write(read.header + '\n')
	outFile.write(read.seq + '\n')
    outFile.close()

def exportPairedFasta(pairedReads, outFileName, outFileNameN):
    outFile = open(outFileName, 'w')
    outFileN = open(outFileNameN, 'w')
    for read in pairedReads:
	outFile.write(read.header1 + '\n')
	outFile.write(read.seq + '\n')
	outFileN.write(read.header1 + '\n')
	outFileN.write(read.seqN + '\n')
    outFile.close()
    outFileN.close()
