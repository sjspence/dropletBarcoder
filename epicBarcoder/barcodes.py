#!/usr/bin/env python

#Within each sample, group by barcode; quantify unique barcode pairings
#Input: Fasta file with droplet barcode, otu, and taxonomic information in the header
#Output: A dictionary where each sampID maps to a dictionary of droplet barcodes:[zOTU1, zOTU2]
#        A dictionary where each zOTU maps to a >80% taxonomy
def createBarcodeDict(inFileName):
    inFile = open(inFileName, 'r')
    barcodeSamples = {}
    taxDict = {}
    for line in inFile:
        if '>' in line:
            line = line.strip().split(';')
            samp = line[0].split('_')[0].replace('>','')
            bc = line[0].split('droplet_bc=')[1]
            otu = line[1]
            tax = line[2].replace('tax=','')
            if samp not in barcodeSamples:
                barcodeSamples[samp] = {bc:[otu]}
            else:
                if bc not in barcodeSamples[samp]:
                    barcodeSamples[samp][bc] = [otu]
                else:
                    barcodeSamples[samp][bc].append(otu)
            if otu not in taxDict:
                taxDict[otu] = tax
    inFile.close()
    return barcodeSamples, taxDict
