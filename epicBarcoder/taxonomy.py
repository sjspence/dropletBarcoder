#!/usr/bin/env python
# 3/27/2017 Sarah J. Spencer, Alm Lab

import pandas as pd

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

#Take in a dictionary of zOTU:taxonomy, then construct a dataframe that
#groups zOTUs by those that share taxonomy (i.e. tOTUs)
#INPUT:  taxonomic dictionary, output of importSintax()
#OUTPUT: pandas dataframe with zOTU indexes, and two columns with unique tOTU
#        ID and full taxonomic string
def tOTUmap(taxDict):
    index = []
    tax_tOTU = {}
    data = []
    i = 1
    for zOTU in taxDict:
        index.append(zOTU)
        tax = taxDict[zOTU]
        if tax not in tax_tOTU:
            tOTU = 'tOtu' + str(i)
            tax_tOTU[tax] = tOTU
            taxList = [tOTU, tax]
            i += 1
        else:
            taxList = [tax_tOTU[tax], tax]
        data.append(taxList)
    columns = ['tOTU', 'taxonomy']
    otuDf = pd.DataFrame(data, index=index, columns=columns)
    return otuDf

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
