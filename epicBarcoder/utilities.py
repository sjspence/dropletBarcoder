import pandas as pd
import os
import subprocess
from collections import defaultdict
from itertools import combinations, chain
from scipy.stats import poisson
from . import io

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


# def move_barcodes_and_type_to_fasta_id(bc_seq, bridge_dict):
#     bridge_dict = {key: ep.expand_primers(ep.reverse_complement(val))
#                    for key, val in bridge_dict.items()}
#     for seq_id, seq in bc_seq:
#         for bridge_id, bridges in bridge_dict.items():
#             for bridge in bridges:
#                 if bridge in seq:
#                     bc, rest = seq.split(bridge)
#                     if len(bc) == 20:
#                         seq_id = "{} barcode={} sequence_type={}".format(seq_id.strip(),
#                                                                          bc, bridge_id)
#                         yield([seq_id, rest])


def process_barcode_info(input_seq_file, output_seq_file, bridge_dict):
    seqs = io.read_fasta(input_seq_file)
    fasta_iter = move_barcodes_and_type_to_fasta_id(seqs, bridge_dict)
    io.write_fasta(fasta_iter, output_seq_file)


def get_len_distr(seqs):
    len_distr_dict = defaultdict(list)
    for seq_id, seq in seqs:
        seq_type_id = seq_id.split("_")[-1]
        len_distr_dict[seq_type_id].append(int(len(seq) / 10) * 10)
    return len_distr_dict


def get_seed_dict(uc_file):
    seed_dict = {}
    with open(uc_file) as f:
        for line in f:
            if line.split("\t")[0] == "H":
                seq_id = line.split("\t")[8].split()[0]
                seed_id = line.split("\t")[9].split()[0]
                seed_dict[seq_id] = seed_id
            if line.split("\t")[0] == "S":
                seq_id = line.split("\t")[8].split()[0]
                seed_dict[seq_id] = seq_id
    return seed_dict


def add_otus_to_fasta(seq_file, output_file, uc_files):
    seeds = {}
    for uc in uc_files:
        seed = get_seed_dict(uc)
        seeds.update(seed)
    seq_acc = []
    for seq_id, seq in io.read_fasta(seq_file):
        short_id = seq_id[1:].split()[0]
        try:
            seed_id = seeds[short_id]
            seq_id = "{} OTU={}".format(seq_id.strip(), seed_id)
        except KeyError:
            pass
        seq_acc.append([seq_id, seq])
    io.write_fasta(seq_acc, output_file)


def fasta_to_bc_otu_table(fasta_file):
    fasta_list = list(io.read_fasta(fasta_file))
    fasta_table = pd.DataFrame([i.strip().split() for i, _ in fasta_list])
    fasta_table['Sample'] = fasta_table[0].str.split("_").apply(lambda x: x[0][1:])
    fasta_table['Barcode'] = fasta_table[5].str.split("=").apply(lambda x: x[1])
    fasta_table['Type'] = fasta_table[6].str.split("=").apply(lambda x: x[1])
    fasta_table[7][fasta_table[7].isnull()] = "OTU=None"
    fasta_table['OTU'] = fasta_table[7].str.split("=").apply(lambda x: x[1])
    fasta_table = fasta_table.loc[:,[0, 'Sample', 'Barcode', 'Type', 'OTU']].rename(columns={0: 'Read'})
    return fasta_table

def process_mapping_file(mapping_file):
    sampIDs = []
    mapping = {}
    readCounts = {}
    with open(mapping_file, 'r') as inFile:
        for line in inFile:
            if '#' not in line:
                line = line.strip().split('\t')
                mapping[line[1]] = line[0].replace('_', 's')
                readCounts[line[1]] = 0
                sampIDs.append(line[0].replace('_', 's'))
    return [sampIDs, mapping, readCounts]


def process_fastq_and_mapping_file(input_file, output_file, mapping_file, quality_summary_file):
    sampIDs, mapping, readCounts = process_mapping_file(mapping_file)
    with open(input_file, 'r') as inFile:
        with open(output_file, 'w') as outFile:
            i = 0
            j = 0
            nextSeq = False
            for line in inFile:
                if nextSeq:
                    outFile.write(line)
                    nextSeq = False
                if i % 4 == 0:
                    for bc in mapping:
                        if bc in line:
                            readCounts[bc] += 1
                            newLine = line.strip().replace('@', '>' + mapping[bc] + '_' + str(j) + ' ')
                            newLine = newLine + ' orig_bc=' + bc + ' new_bc=' + bc + ' bc_diffs=0\n'
                            outFile.write(newLine)
                            nextSeq = True
                            j += 1
                i += 1
    total = 0
    with open(quality_summary_file, 'w') as summaryFile:
        for s in sampIDs:
            for bc in mapping:
                if mapping[bc] == s:
                    summaryFile.write(s + '\t' + str(readCounts[bc]) + '\n')
                    total += readCounts[bc]
        summaryFile.write('Total\t' + str(total))


def filter_significant_connections(connection_file, abundance_file, sig_above_file, sig_below_file):
    conn = pd.read_csv(connection_file, header=None)
    abu = pd.read_csv(abundance_file, header=None)
    otu_count = abu[1].sum()
    abu[2] = abu[1]/otu_count
    pairwise = pd.DataFrame(list(combinations(list(abu[0]), 2)))
    pairwise_abu1 = pd.merge(pairwise, abu, on=0)
    pairwise_abu2 = pd.merge(pairwise_abu1, abu, left_on='1_x', right_on=0)
    pairwise_abu = pairwise_abu2.loc[:, ['0_x', '1_x', '2_x', '2_y']]
    pairwise_abu.columns = ['Left', 'Right', 'Left_perc', 'Right_perc']
    pairwise_abu['Lambda'] = pairwise_abu['Left_perc'] * pairwise_abu['Right_perc'] * otu_count
    pairwise_tot = pd.merge(pairwise_abu, conn, left_on=['Left', 'Right'], right_on=[0, 1])
    pairwise_tot = pairwise_tot.loc[:, ['Left', 'Right', 'Lambda', 2]]
    pairwise_tot.columns = ['Left', 'Right', 'Lambda', 'Observed']
    pairwise_tot['p-val'] = poisson.pmf(pairwise_tot['Observed'], pairwise_tot['Lambda'])
    pairwise_sig_above = pairwise_tot[(pairwise_tot['p-val'] < 0.001) & (pairwise_tot['Observed'] > pairwise_tot['Lambda'])]
    pw_s_a_out = pairwise_sig_above.loc[:, ['Left', 'Right', 'Observed']]
    pw_s_a_out['Color'] = '#ff0000'
    pw_s_a_out['Label'] = 'Label'
    pw_s_a_out.to_csv(sig_above_file, index=None) 
    pairwise_sig_below = pairwise_tot[(pairwise_tot['p-val'] < 0.001) & (pairwise_tot['Observed'] < pairwise_tot['Lambda'])]
    pw_s_b_out = pairwise_sig_below.loc[:, ['Left', 'Right', 'Observed']]
    pw_s_b_out['Color'] = '#000099'
    pw_s_b_out['Label'] = 'Label'
    pw_s_b_out.to_csv(sig_below_file, index=None) 


def output_abunds_and_connections(input_file, abund_output_file, connection_output_file):
    a = pd.read_csv(input_file, sep="\t", header=None)
    a[1] = a[1].str.split(",")
    singlets = pd.Series([i[0] for i in list(a[~a[1].isnull()][1]) if len(i) == 1])
    singlets.value_counts().reset_index().to_csv(abund_output_file, index=None, header=None)
    combs = list(pd.Series([list(combinations(i, 2)) for i in list(a[~a[1].isnull()][1]) if len(i) > 1]))
    connections = pd.DataFrame([sorted(i) for i in list(chain.from_iterable(combs))])
    connections = connections.groupby(0)[1].value_counts()
    connections.name = 'Count'
    connections = connections.reset_index()
    connections['Color'] = '#ff0000'
    connections['Label'] = 'Label'
    connections.to_csv(connection_output_file, header=None, index=None)


def output_functions(input_file, output_file, gene_name):
    a = pd.read_csv(input_file, sep="\t", header=None)
    b = a[~(a[1].isnull()) & ~(a[3].isnull()) & ~(a[1].str.contains(",", na=False)) & ~(a[3].str.contains(",", na=False))]
    b.columns = ['Barcode', 'OTU', 'Euk_OTU', 'Func']
    c = b.groupby('OTU')['Func'].value_counts()
    c.name = 'Count'
    c = c.reset_index()
    d = pd.pivot_table(c, values='Count', index='OTU', columns='Func', fill_value=-1)
    d[d>1] = 1
    d[gene_name].to_csv(output_file, header=None)
