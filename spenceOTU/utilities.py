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
#    original_fasta_list = [[seq_id.strip().split(">")[1], seq]
#                           for seq_id, seq
#                           in read_fasta(seq_name)]
#    uc_table = []
#    with open(seq_name + ".uc") as f:
#        for line in f:
#            if line[0] != "#":
#                uc_table.append(line.strip().split("\t"))
#    clustered_fasta_list = [[seq_id.strip().split(">")[1], seq]
#                            for seq_id, seq
#                            in read_fasta(clustered_output)]
#    os.remove(seq_name + ".uc")
#    os.remove(clustered_output)
#    return uc_table, clustered_fasta_list, original_fasta_list
