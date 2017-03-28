#!/usr/bin/env python
# 3/27/2017 Sarah J. Spencer, Alm Lab

#Import SINTAX taxonomy
#Input: usearch SINTAX taxonomy file name
#Output: python dictionary mapping otu ID to 
#        list [taxonomy probabilities, final taxonomy]
def importSintax(inFileName):
    taxDict = {}
    inFile = open(inFileName, 'r')
    for line in inFile:
        line = line.strip().split('\t')
        otuID = line[0].split(' ')[0]
        taxProbs = line[1]
        tax = line[3]
        taxDict[otuID] = [taxProbs, tax]
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
