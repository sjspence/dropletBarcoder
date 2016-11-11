import os
import subprocess

def getExtension(fileName):
    fileName = fileName.split('.')
    fileExt = '.' + fileName[len(fileName)-1]
    return fileExt

def clusterWithUsearch(usearchPath, inFile, percIdentity):
    env = os.environ
    percIdentity = str(percIdentity)
    extension = getExtension(inFile)
    outFile = inFile.replace(extension, '_' + percIdentity + '.fa')
    clustFile = inFile.replace(extension, '_' + percIdentity + '.uc')
    subprocess.call([usearchPath, "-cluster_fast", inFile,
                     "-sort", "length", "-id", percIdentity,
                     "-centroids", outFile, "-uc",
                     clustFile], env=env)
