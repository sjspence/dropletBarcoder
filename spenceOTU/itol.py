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
