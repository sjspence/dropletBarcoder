#!/usr/bin/env python
# 3/27/2017 Sarah J. Spencer, Alm Lab

#Import usearch -cluster_fast hits
#Input: usearch '.uc' file
#Output: python dictionary mapping seeds: [list of hits]
def importClusterFast(usearchFileName):
    usearchFile = open(usearchFileName, 'r')
    seedDict = {}
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
