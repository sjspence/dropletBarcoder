#!/usr/bin/env python
# 3/27/2017 Sarah J. Spencer, Alm Lab

#Import SINTAX taxonomy
#INPUT:  usearch SINTAX taxonomy file name
#	 'full' or 'final', see below
#OUTPUT: python dictionary mapping otu ID to 
#        'full': list [taxonomy probabilities, final taxonomy]
#	 'final': string with >80% taxonomic assignments
def importSintax(inFileName, completeness):
    taxDict = {}
    inFile = open(inFileName, 'r')
    for line in inFile:
        line = line.strip().split('\t')
        otuID = line[0].split(' ')[0].split(';')[0]
        taxProbs = line[1]
        tax = line[3]
	if completeness == 'full':
            taxDict[otuID] = [taxProbs, tax]
	elif completeness == 'final':
	    taxDict[otuID] = tax
    inFile.close()
    return taxDict

#Import mothur taxonomy file
#Input: taxonomy file name
#	taxType = 'mothur' or ... could add other functionality
#Output: 'mothur': python dictionary mapping read/otu ID to [list of 
#	 taxonomies separated by semicolon, probability of correct assignment]
def importTaxonomy(taxonomyFileName, taxType):
    taxFile = open(taxonomyFileName, 'r')
    taxDict = {}
    for line in taxFile:
        line = line.strip().split('\t')
        if taxType == 'mothur':
            taxDict[line[0]] = line[1:]
    taxFile.close()
    return taxDict
