#!/usr/bin/env python

import os

def fixKey(entry):
    if "_" in entry:
	#entry = entry.replace("_"," ")
        entry = "\'" + entry + "\'"
    return entry

def itolHeatmap(df, outFileName):
    outFile = open(outFileName, 'w')
    outFile.write('DATASET_HEATMAP\n')
    outFile.write('SEPARATOR COMMA\n')
    outFile.write('DATASET_LABEL,abundances heat\n')
    outFile.write('COLOR,#ff0000\n')
    fields = ','.join(df.columns.values.tolist())
    outFile.write('FIELD_LABELS,' + fields + '\n')
    outFile.write('COLOR_MIN,#ff0000\n')
    outFile.write('COLOR_MAX,#0000ff\n')
    outFile.write('DATA\n')
    for index, row in df.iterrows():
        outFile.write(index + ',')
        abundList = []
        for samp in df.columns.values.tolist():
            abundList.append(str(row[samp]))
        outFile.write(','.join(abundList) + '\n')
    outFile.close()

#Create iTol simple bar files for each column in a dataframe
#INPUT:  dataframe, with e.g. sample IDs in headers
#	 directory to store all the files generated
#OUTPUT: files written to outDirectory, outDirectory also zipped in the same
#	     location
def itolSimpleBar(df, outDirectory):
    if not os.path.exists(outDirectory):
        os.makedirs(outDirectory)
    if outDirectory[-1] != '/':
        outDirectory = outDirectory + '/'
    for samp in df:
        outFile = open(outDirectory + samp + '_abund.txt', 'w')
        outFile.write('DATASET_SIMPLEBAR\nSEPARATOR COMMA\n')
        outFile.write('DATASET_LABEL,' + samp + '_abund\n')
        outFile.write('COLOR,#547261\n')
	outFile.write('WIDTH,100\n')
        outFile.write('DATA\n')
        for otu in list(df.index.values):
            outFile.write(otu + ',' + str(df[samp][otu]) + '\n')
        outFile.close()
    outZip = outDirectory[0:len(outDirectory)-1] + '.zip'
    os.system('zip -r ' + outZip + ' ' + outDirectory)

#Create iTol connection files for each column in a dataframe, where the indexes
#are otu pairs separated by a double underscore
#INPUT:  dataframe with samples as columns, otu pairs separated by double
#            underscore as indexes, and number of pairs observed as data
#	 outDirectory: directory to write files to, then zip
#	 tag: a string that describes the data (e.g. 'all', 'significant')
#	 color: hex string used for connection line colors
#OUTPUT: files written to outDirectory, outDirectory also zipped in the same
#            location
def itolConnections(pairDf, outDirectory, tag, color):
    if not os.path.exists(outDirectory):
        os.makedirs(outDirectory)
    if outDirectory[-1] != '/':
        outDirectory = outDirectory + '/'
    for samp in pairDf:
        outFile = open(outDirectory + samp + '_' + tag + '.txt', 'w')
        outFile.write('DATASET_CONNECTION\n')
        outFile.write('SEPARATOR COMMA\n')
        outFile.write('DATASET_LABEL,' + samp + '_' + tag + '\n')
        outFile.write('COLOR,' + color + '\n')
        outFile.write('DRAW_ARROWS,0\n')
        outFile.write('MAXIMUM_LINE_WIDTH,10\n')
        outFile.write('CENTER_CURVES,1\nALIGN_TO_LABELS,1\n')
        outFile.write('DATA\n')
        for pair in list(pairDf.index.values):
            if pairDf[samp][pair] == 0.0:
                continue
            otu1 = pair.split('__')[0]
            otu2 = pair.split('__')[1]
            outFile.write(otu1 + ',' + otu2 + ',' + str(pairDf[samp][pair]) + \
			  ',' + color + '\n')
        outFile.close()
    outZip = outDirectory[0:len(outDirectory)-1] + '.zip'
    os.system('zip -r ' + outZip + ' ' + outDirectory)
