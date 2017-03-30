#!/usr/bin/env python
# 10/6/2016 Sarah J. Spencer, Alm Lab

import pandas as pd

#Take in a dictionary of zOTU:taxonomy, then construct a dataframe that
#groups zOTUs by those that share taxonomy (i.e. tOTUs)
#INPUT:  taxonomic dictionary, output of importSintax()
#OUTPUT: pandas dataframe with zOTU indexes, and two columns with unique tOTU
#	 ID and full taxonomic string
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
