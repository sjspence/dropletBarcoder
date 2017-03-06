import pandas as pd
import os
import subprocess
from collections import defaultdict
from itertools import combinations, chain
from scipy.stats import poisson
from . import io
from . import dereplicate
from . import reads

env = os.environ

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
            type_id = seq_id.strip().split()[-1].split("=")[1]
            seq_id = "{} OTU={}__{}".format(seq_id.strip(), type_id, seed_id)
        except KeyError:
            seed_id = "unclassified"
            type_id = seq_id.strip().split()[-1].split("=")[1]
            seq_id = "{} OTU={}__{}".format(seq_id.strip(), type_id, seed_id)
        seq_acc.append([seq_id, seq])
    io.write_fasta(seq_acc, output_file)


def fasta_to_bc_otu_table(fasta_file, output_file=None):
    fasta_list = list(io.read_fasta(fasta_file))
    fasta_table = pd.DataFrame([i.strip().split() for i, _ in fasta_list])
    fasta_table['Sample'] = fasta_table[0].str.split("_").apply(lambda x: x[0][1:])
    fasta_table['Barcode'] = fasta_table[5].str.split("=").apply(lambda x: x[1])
    fasta_table['Type'] = fasta_table[6].str.split("=").apply(lambda x: x[1])
    fasta_table[7][fasta_table[7].isnull()] = "OTU=None"
    fasta_table['OTU'] = fasta_table[7].str.split("=").apply(lambda x: x[1])
    fasta_table = fasta_table.loc[:, [0, 'Sample', 'Barcode',
                                      'Type', 'OTU']].rename(columns={0: 'Read'})
    if output_file:
        fasta_table.to_csv(output_file, index=None)
    return fasta_table


def get_grouped_table(fasta_table):
    if isinstance(fasta_table, str):
        fasta_table = pd.read_csv(fasta_table)
    grouped_table = fasta_table.groupby(['Sample', 'Barcode'])['OTU'].apply(set)
    filtered_table = grouped_table.apply(lambda x: [i for i in list(x) if (i != "16S__unclassified")
                                                    and (i != "18S__unclassified")])
    filtered_table = filtered_table[filtered_table.apply(lambda x: len(x) > 0)]
    return filtered_table


def get_singletons_and_connections(grouped_table):
    bacterial_otus = grouped_table.apply(lambda x: [i for i in x if "16S__" in i])
    bacterial_otus = bacterial_otus[bacterial_otus.apply(lambda x: len(x) > 0)]

    singletons = bacterial_otus[bacterial_otus.apply(lambda x: len([i for i in x if "16S__" in i]) == 1)]\
                 .apply(lambda x: x[0])
    singletons = singletons.groupby(level='Sample').value_counts()
    singletons.name = 'Count'
    singletons = singletons.reset_index()
    singletons['OTU'] = singletons['OTU'].str.split("__").apply(lambda x: x[1])

    multiples = bacterial_otus[bacterial_otus.apply(lambda x: len([i for i in x if "16S__" in i]) > 1)]
    acc = []
    for entry_id, entry in multiples.groupby(level='Sample'):
        tmp = [combinations(i, 2) for i in list(entry)]
        tmp = [sorted(i) for i in chain.from_iterable(tmp)]
        tmp = pd.DataFrame(tmp)
        tmp = tmp.applymap(lambda x: x.split("__")[1])
        tmp = tmp.groupby(0)[1].value_counts()
        tmp.name = 'Count'
        tmp = tmp.reset_index()
        tmp['Sample'] = entry_id
        tmp = tmp.set_index('Sample')
        acc.append(tmp)
    multiples = pd.concat(acc)
    multiples = multiples.reset_index()

    return [singletons, multiples]


def write_connections_and_abundances(conn_abund_list, file_name_prefix = None):
    abundances, connections = conn_abund_list
    for sample, obs in abundances.groupby('Sample'):
        sample_name = sample + "_abunds.csv"
        if file_name_prefix:
            sample_name = file_name_prefix + sample_name
        obs.iloc[:, [1,2]].to_csv(sample_name, index=None)

    for sample, obs in connections.groupby('Sample'):
        sample_name = sample + "_connections.csv"
        if file_name_prefix:
            sample_name = file_name_prefix + sample_name
        obs.to_csv(sample_name, index=None)


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


def make_otus_and_assign(input_file, db_dir, usearchPath):
    noPrimerReads = reads.importFasta(input_file)
    uniqueDict = dereplicate.getUniqueSeqs(noPrimerReads, '05_unique_seqs.fasta')
    subprocess.call([usearchPath, '-unoise2', '05_unique_seqs.fasta', '-fastaout', '06_denoised.fa',
                    '-otudbout', '06_db.fa', '-minampsize', '1'], env=env)
    outFile = open(db_dir + 'HOMD_16S_rRNA_RefSeq_V14.51.p9_sintax.fasta', 'w')
    taxDict = {}
    with open(db_dir + 'HOMD_16S_rRNA_RefSeq_V14.51.qiime.taxonomy', 'r') as t:
        for line in t:
            line = line.strip().split('\t')
            taxID = line[0]
            tax = line[1].strip().replace('__',':')
            tax = tax.replace(';',',')
            taxDict[taxID] = tax
    with open(db_dir + 'HOMD_16S_rRNA_RefSeq_V14.51.p9.fasta', 'r') as f:
        for line in f:
            if '>' in line:
                line = line.strip().split(' ')
                taxInfo = taxDict[line[0].replace('>','')]
                outLine = line[0] + ';tax=' + taxInfo + ';'
                for i in line:
                    if 'HOT' in i:
                        outLine += i + ';'
                outFile.write(outLine + '\n')
            else:
                outFile.write(line)
    outFile.close()
    subprocess.call([usearchPath, '-makeudb_sintax', db_dir + 'HOMD_16S_rRNA_RefSeq_V14.51.p9_sintax.fasta',
                    '-output', db_dir + 'HOMD_16S_rRNA_RefSeq_V14.51.p9_sintax.udb'], env=env)
    subprocess.call([usearchPath, '-sintax', '06_denoised.fa',
                    '-db', db_dir + 'HOMD_16S_rRNA_RefSeq_V14.51.p9_sintax.udb',
                    '-tabbedout', '07_denoised.sintax',
                    '-strand', 'plus', '-sintax_cutoff', '0.8', '-threads', '4'], env=env)
    denoised = reads.importFasta('06_denoised.fa')
    taxDict = io.importSintax('07_denoised.sintax')
    dereplicate.otuToHeaders(denoised, taxDict, uniqueDict, '08_all_seqs_tax.fa')


def filter_significant_connections(connection_file, abundance_file, sig_above_file, sig_below_file):
    conn = pd.read_csv(connection_file)
    conn = conn.iloc[:, 1:]
    conn.columns = [0, 1, 2]
    abu = pd.read_csv(abundance_file)
    abu = abu.iloc[:, 1:]
    abu.columns = [0, 1]
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


def process_unoise_fasta(input_file, output_file):
    fst = io.read_fasta(input_file)
    acc = []
    for seq_id, seq in fst:
        seq_id_split = seq_id.split()
        bc_field = seq_id_split[-1]
        bc, otu = bc_field.split(";")[:2]
        bc = "barcode=" + bc.split("=")[1]
        new_seq_id = " ".join(seq_id_split[:-1]) + " " + bc + " sequence_type=16S " + "OTU=" + otu
        acc.append([new_seq_id, seq])
    io.write_fasta(acc, output_file)
