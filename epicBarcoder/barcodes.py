#!/usr/bin/env python

#Within each sample, group by barcode; assign each barcode a list of zOTUs
#for each read that contains that barcode (i.e. list can contain redundancies)
#INPUT:  Fasta file with droplet barcode and zOTU information in the header
#OUTPUT: A dictionary where each sampID maps to a dictionary of droplet 
#	 barcodes:[zOTU1, zOTU1, zOTU2, ...]
def createBarcodeDict(inFileName):
    inFile = open(inFileName, 'r')
    barcodeSamples = {}
    for line in inFile:
        if '>' in line:
            line = line.strip().split(';')
            samp = line[0].split('_')[0].replace('>','')
            bc = line[0].split('droplet_bc=')[1]
            otu = line[1]
            if samp not in barcodeSamples:
                barcodeSamples[samp] = {bc:[otu]}
            else:
                if bc not in barcodeSamples[samp]:
                    barcodeSamples[samp][bc] = [otu]
                else:
                    barcodeSamples[samp][bc].append(otu)
    inFile.close()
    return barcodeSamples
