
def makeTreeConstraint(inFileName, outFileName):
    inFile = open(inFileName, 'r')
    taxAssignments = []
    for line in inFile:
        x = line.strip().split('\t')
        taxAssignments.append(x)
    inFile.close()
    seqID = []
    tax = []
    for i in taxAssignments:
        seqID.append(i[0])
        tax.append(i[1].split(';'))
    taxSet = [[] for i in range(7)]    #kingdom, phylum, class, order, family, genus, species
    for a in tax:
        for i in range(7):
            if len(a) > i and a[i] not in taxSet[i]:
                taxSet[i].append(a[i])
    #now make a list of each binary identifier
    outFile = open(outFileName, 'w')
    test = []
    posNum = []
    for x in taxAssignments:
        taxBinary = [[] for i in range(7)]
        pos = taxAssignments.index(x)
        test.append(seqID[pos])
        posNum.append(pos)
        for i in range(7):
            for j in taxSet[i]:
                k = "0"
                if j in x[1]:
                    k = "1"
                taxBinary[i].append(k)
        tBins=[]
        for w in taxBinary:
            tBins.append(''.join(w))
        outFile.write(">" + seqID[pos] + "\n" + ''.join(tBins) + "\n")
    outFile.close()
