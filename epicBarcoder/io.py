
from itertools import groupby
import pandas as pd

import pandas as pd

def exportPairedFasta(pairedReads, outFileName, outFileNameN):
    outFile = open(outFileName, 'w')
    outFileN = open(outFileNameN, 'w')
    for read in pairedReads:
        outFile.write(read.header1 + '\n')
        outFile.write(read.seq + '\n')
        outFileN.write(read.header1 + '\n')
        outFileN.write(read.seqN + '\n')
    outFile.close()
    outFileN.close()

def exportOTUtable(df, taxDict, taxType, outFileName):
    outFile = open(outFileName, 'w')
    outFile.write('\t')
    outFile.write('\t'.join(df.columns.values))
    if taxType == 'mothur':
        outFile.write('\ttaxonomy\tprobability\n')
    else:
        outFile.write('\ttaxonomy\n')
    for row in df.iterrows():
        index, data = row
        outFile.write(index + '\t')
        outFile.write('\t'.join(map(str, data.tolist())) + '\t')
        outFile.write('\t'.join(taxDict[index]) + '\n')

#Include the path to a text file output from a Qiime OTU table and
#a true/false statement indicating if there is taxonomic information
def importQiimeOTU(textFileName, taxIncluded):
    textFile = open(textFileName, 'r')
    i = 0 #first line of file is a note, second line is header data
    header = []
    data = []
    index = []
    for line in textFile:
        lineList = line.strip().split('\t')
        if i == 1:
            lineList = lineList[1:] #remove otu designation from first position
            header = lineList
        elif i > 1:
            index.append(lineList[0])
            lineList = lineList[1:]
            if taxIncluded:
                data.append(map(float, lineList[0:-1]) + [lineList[-1]])
            else:
                data.append(map(float, lineList))
        i += 1
    textFile.close()
    df = pd.DataFrame(data, index=index, columns=header)
    return df


def read_fasta(fasta):
    """ Read a fasta file and yield pairs of fasta ids and DNA sequences.
    """
    with open(fasta) as f:
        grouped = groupby(f, lambda x: x[0] == ">")
        for cond, entry in grouped:
            if cond:
                fasta_id = next(entry)
                _, seq_iter = next(grouped)
                seq = ''.join([line.strip() for line in seq_iter]).upper()
                yield([fasta_id, seq])

def write_fasta(fasta_iter, output_file, size_limit=100000):
    """ Write fasta iterables (such as those created by create_fasta_iter) into fasta files.
    """
    with open(output_file, "w") as f:
        for fasta_id, seq in fasta_iter:
            if fasta_id[0] != ">":
                fasta_id = ">" + fasta_id
            if "\n" not in fasta_id:
                fasta_id += "\n"
            if len(seq) < size_limit:
                output = "{}{}\n".format(fasta_id, seq)
                f.write(output)


def read_fastq(fastq):
    """
    Read a fastq file and yield chunks of four lines, corresponding to each entry
    :param fastq: str; fastq file name
    :return: list; a list of id, sequence, other id and quality score
    """
    with open(fastq) as f:
        for line in f:
            id_line = line
            seq_line = next(f)
            other_id_line = next(f)
            quality_score = next(f)
            yield([id_line, seq_line, other_id_line, quality_score])


def write_fastq(fastq_iter, output_file, size_limit=100000):
    """ Write fastq iterables (such as those created by create_fastq_iter) into fastq files.
    """
    with open(output_file, "w") as f:
        for fastq_id, seq, other_id, quality in fastq_iter:
            if len(seq) < size_limit:
                fastq_id = fastq_id.strip()
                seq = seq.strip()
                other_id = other_id.strip()
                quality = quality.strip()
                out_str = "{}\n{}\n{}\n{}\n".format(fastq_id, seq, other_id, quality)
                f.write(out_str)
