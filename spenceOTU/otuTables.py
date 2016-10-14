#!/usr/bin/env python
# 10/6/2016 Sarah J. Spencer, Alm Lab

import pandas as pd

def importUsearchHits(usearchFileName):
    usearchFile = open(usearchFileName, 'r')
    seedDict = {} #seeds: [list of hits]
    for line in usearchFile:
	line = line.strip().split('\t')
	if line[0] == 'S':
	    if line[8] not in seedDict:
		seedDict[line[8]] = []
	elif line[0] == 'H':
	    if line[9] not in seedDict:
		seedDict[line[9]] = [line[8]]
	    else:
		seedDict[line[9]].append(line[8])
    usearchFile.close()
    return seedDict

def importTaxonomy(taxonomyFileName, taxType):
    taxFile = open(taxonomyFileName, 'r')
    taxDict = {}
    for line in taxFile:
	line = line.strip().split('\t')
	if taxType == 'mothur':
	    taxDict[line[0]] = line[1:]
    return taxDict
    taxFile.close()

def buildOTUtable(sampleIDs, usearchHits):
    seeds = []
    otuTable = []
    for seed in usearchHits:
	seeds.append(seed)
	outLine = []
	for samp in sampleIDs:
	    ct = 0
	    if samp in seed:
		ct += 1
	    for read in usearchHits[seed]:
		if samp in read:
		    ct += 1
	    outLine.append(ct)
	otuTable.append(outLine)
    otuPanda = pd.DataFrame(otuTable, index=seeds, columns=sampleIDs)
    return otuPanda

#Get a dictionary of read:seed, matching read seq_names with their OTU seed sequence
def invertHits(hits):
    readToOTU = {}
    for h in hits:
        readToOTU[h] = h
        for seqID in hits[h]:
            readToOTU[seqID] = h
    return readToOTU
