#!/usr/bin/env python
# 10/6/2016 Sarah J. Spencer, Alm Lab

import pandas as pd

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
