import os
import subprocess
import random
import string
import pickle
from collections import defaultdict
from itertools import combinations
from scipy.stats import poisson
import pandas as pd
from . import io
from . import dereplicate
from . import reads


env = os.environ


connection_template = '''DATASET_CONNECTION
SEPARATOR COMMA
DATASET_LABEL,{}
COLOR,#ff0ff0
DRAW_ARROWS,0
ARROW_SIZE,0
MAXIMUM_LINE_WIDTH,10
CURVE_ANGLE,0
CENTER_CURVES,1
ALIGN_TO_LABELS,0
DATA
#NODE1,NODE2,WIDTH,COLOR,LABEL
'''


abund_hist_template = '''DATASET_SIMPLEBAR
SEPARATOR COMMA
DATASET_LABEL,{}
WIDTH,300
BORDER_WIDTH,5
COLOR,{}
DATA
#ID1,value1
'''


popup_template = '''POPUP_INFO
SEPARATOR COMMA
DATA
#NODE_ID,POPUP_TITLE,POPUP_CONTENT
#popup for leaf node 9606
#9606,Homo sapiens info popup,Detailed stuff
'''


def generate_id(size=6):
    """ Generate random sequences of characters for temporary file names.
    """
    chars = string.ascii_uppercase + string.digits
    return ''.join(random.choice(chars) for _ in range(size))


def getExtension(fileName):
    fileName = fileName.split('.')
    fileExt = '.' + fileName[len(fileName)-1]
    return fileExt


def clusterWithUsearch(usearchPath, inFile, percIdentity):
    percIdentity = str(percIdentity)
    extension = getExtension(inFile)
    outFile = inFile.replace(extension, '_' + percIdentity + '.fa')
    clustFile = inFile.replace(extension, '_' + percIdentity + '.uc')
    subprocess.call([usearchPath, "-cluster_fast", inFile,
                     "-sort", "length", "-id", percIdentity,
                     "-centroids", outFile, "-uc",
                     clustFile], env=env)


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


def parse_unoise(unoise_file, seq_type):
    fasta_iter = io.read_fasta(unoise_file)
    otu_acc = []
    tax_acc = []
    for line, _ in fasta_iter:
        sample = line.split()[0].split("_")[0][1:]
        bc_and_otu = line.split()[5].split("=")[1]
        tax = line.split()[5].split("=")[2].split(",")
        bc = bc_and_otu.split(";")[0]
        otu = bc_and_otu.split(";")[1]
        otu_acc.append([sample, bc, seq_type, otu])
        tax_acc.append([otu] + tax)
    otu = pd.DataFrame(otu_acc)
    otu.columns = ['Sample', 'Barcode', 'Type', 'OTU']
    grouped_table = otu.groupby(['Sample', 'Barcode', 'Type'])['OTU'].apply(list)
    tax = pd.DataFrame(tax_acc)
    tax.columns = ['OTU', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    return [grouped_table, tax]


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
    grouped_table = fasta_table.groupby(['Sample', 'Barcode', 'Type'])['OTU'].apply(list)
    return [grouped_table, []]


def get_singletons(grouped_table):
    singletons = grouped_table[grouped_table.apply(lambda x: len(x) == 1)]
    singletons = singletons.reset_index()
    singletons['OTU'] = singletons['OTU'].apply(lambda x: x[0])
    singletons = singletons.loc[:, ['Sample', 'OTU']]
    return singletons


def get_connections(grouped_table):
    print("Filtering out singletons..")
    connections = grouped_table[grouped_table.apply(lambda x: len(set(x)) > 1)]
    connections = connections.reset_index()
    print("Filtering out multiple entries..")
    connections['OTU'] = connections['OTU'].apply(lambda x: sorted(list(set(x))))
    return connections


def expand_connections(connections):
    acc = []
    connections['OTU'] = connections['OTU'].apply(lambda x: combinations(x, 2))
    for row in connections.iterrows():
        sample = row[1]['Sample']
        bc = row[1]['Barcode']
        otu_type = row[1]['Type']
        new_frame = pd.DataFrame([sorted(j) for j in row[1]['OTU']])
        new_frame['Sample'] = sample
        new_frame['Barcode'] = bc
        new_frame['Type'] = otu_type
        acc.append(new_frame)
    expanded = pd.concat(acc)
    expanded['Connection'] = expanded[0] + "," + expanded[1]
    expanded = expanded.loc[:, ['Sample', 'Connection']]
    expanded.index = range(len(expanded.index))
    return expanded


def grouper(input_fasta, unoise, seq_type=None):
    print("Reading in fasta file..")
    if unoise:
        tables = parse_unoise(input_fasta, seq_type)
    else:
        tables = get_grouped_table(fasta_to_bc_otu_table(input_fasta))
    return tables

def filter_significant_connections(conn, abu):
    conn.columns = [0, 1, 2]
    abu.columns = [0, 1]
    otu_count = abu[1].sum()
    abu[2] = abu[1]/otu_count
    pairwise = pd.DataFrame([sorted(i) for i in combinations(list(abu[0]), 2)])
    pairwise_abu1 = pd.merge(pairwise, abu, on=0)
    pairwise_abu2 = pd.merge(pairwise_abu1, abu, left_on='1_x', right_on=0)
    pairwise_abu = pairwise_abu2.loc[:, ['0_x', '1_x', '2_x', '2_y']]
    pairwise_abu.columns = ['Left', 'Right', 'Left_perc', 'Right_perc']
    pairwise_abu['Lambda'] = pairwise_abu['Left_perc'] * pairwise_abu['Right_perc'] * otu_count
    pairwise_tot = pd.merge(pairwise_abu, conn, left_on=['Left', 'Right'], right_on=[0, 1])
    pairwise_tot = pairwise_tot.loc[:, ['Left', 'Right', 'Lambda', 2]]
    pairwise_tot.columns = ['Left', 'Right', 'Lambda', 'Observed']
    pairwise_tot['p-val'] = poisson.pmf(pairwise_tot['Observed'], pairwise_tot['Lambda'])
    return pairwise_tot


class BarcodeContainer(object):

    def __init__(self, input_16S=None, input_18S=None, input_funcs=None, unoise=False):
        ''' Input after preparing the fasta ids using function add_otus_to_fasta'''
        self.type_dict = {}
        if input_16S:
            print("Parsing 16S data..")
            self.type_dict['16S'], self.tax = grouper(input_16S, unoise, '16S')
            self.bact_connections = self.__get_connections('16S')
            self.bact_singletons = self.__get_singletons('16S')
        if input_18S:
            print("Parsing 18S data..")
            self.type_dict['18S'], _ = grouper(input_18S, unoise, '18S')
            self.euk_connections = self.__get_connections('18S')
            self.euk_singletons = self.__get_singletons('18S')
        if input_funcs:
            print("Parsing functional gene data..")
            self.type_dict['Funcs'] = grouper(input_funcs, unoise, 'Func')

        print("Extracting sample names..")
        self.samples = self.__get_samples()

        print("Pickling self..")
        pickle_name = input_16S.split(".")[0] + ".pickle"
        with open(pickle_name, "wb") as f:
            pickle.dump(self, f)

        print("Writing out bacterial iTOL files..")
        self.write_itol_files('16S')

    def __get_singletons(self, seq_type):
        print("Extracting singletons..")
        table = self.type_dict[seq_type]
        singletons = get_singletons(table)
        return singletons

    def __get_connections(self, seq_type):
        table = self.type_dict[seq_type]
        connections = get_connections(table)
        print("Expanding connection combinations..")
        expanded = expand_connections(connections)
        return expanded

    def __get_samples(self):
        samples = pd.concat([self.type_dict[i] for i in ['16S', '18S', 'Funcs']
                             if i in self.type_dict])
        samples = samples.reset_index()['Sample']
        samples = sorted(set(samples))
        return samples

    def write_itol_files(self, seq_type):
        popup_name = samples[0] + "_popup.txt"
        self.get_popup_for_itol(out_file=popup_name)
        for sample in self.samples:
            file_type = 'abunds'
            file_name = "{}_{}_{}.txt".format(sample, seq_type, file_type)
            print("Writing abundance file {}".format(file_name))
            self.get_itol_abunds(seq_type, sample, color="#ff0000", out_file=file_name, label=sample + "_abund")

            file_type = 'tot_connections'
            file_name = "{}_{}_{}.txt".format(sample, seq_type, file_type)
            print("Writing total connection file {}".format(file_name))
            self.get_total_itol_connections(seq_type, sample, color="#9aa0a6", out_file=file_name, label=sample + "_tot")

            file_type = 'below_connections'
            file_name = "{}_{}_{}.txt".format(sample, seq_type, file_type)
            print("Writing lower connection file {}".format(file_name))
            self.get_itol_sig_below_connections(seq_type, sample, color="#0000ff", out_file=file_name, label=sample + "_below")

            file_type = 'above_connections'
            file_name = "{}_{}_{}.txt".format(sample, seq_type, file_type)
            print("Writing elevated connection file {}".format(file_name))
            self.get_itol_sig_above_connections(seq_type, sample, color="#ff0000", out_file=file_name, label=sample + "_above")

    def get_total_itol_connections(self, seq_type, sample, out_file=None,
                                   color='#9AA0A6', label='Label'):
        if seq_type == '16S':
            conns = self.bact_connections
        elif seq_type == '18S':
            conns = self.euk_connections
        else:
            raise TypeError('Unknown sequence type')
        conns = conns[conns['Sample'] == sample]
        conns = conns['Connection'].value_counts().reset_index()
        conns['Connection'] = conns['Connection'].apply(str)
        conns['Color'] = color
        conns['Label'] = label
        conns = conns['index'] + "," + conns['Connection'] + "," \
                + conns['Color'] + "," + conns['Label']
        conns = "\n".join(sorted(list(conns)))
        conn_template = connection_template.format(label)
        conns = conn_template + conns
        if out_file:
            with open(out_file, 'w') as f:
                f.write(conns)
        else:
            return conns

    def get_itol_abunds(self, seq_type, sample, color, out_file=None, label='Label'):

        if seq_type == '16S':
            singletons = self.bact_singletons
        elif seq_type == '18S':
            singletons = self.euk_singletons
        else:
            raise TypeError('Unknown sequence type')
        singletons = singletons[singletons['Sample'] == sample]
        singletons = singletons['OTU'].value_counts().reset_index()
        singletons['OTU'] = singletons['OTU'].apply(str)
        singletons = sorted(list(singletons['index'] + "," + singletons['OTU']))
        singletons = "\n".join(singletons)
        hist_template = abund_hist_template.format(label, color)
        singletons = hist_template + singletons
        if out_file:
            with open(out_file, 'w') as f:
                f.write(singletons)
        else:
            return singletons

    def get_significant_connections(self, seq_type, sample):
        if seq_type == '16S':
            singletons = self.bact_singletons
            conns = self.bact_connections
        if seq_type == '18S':
            singletons = self.euk_singletons
            conns = self.euk_connections
        singletons = singletons[singletons['Sample'] == sample]
        singletons = singletons['OTU'].value_counts().reset_index()
        conns = conns[conns['Sample'] == sample]
        conns = conns['Connection'].value_counts().reset_index()
        conns['Left'] = conns['index'].str.split(",").apply(lambda x: x[0])
        conns['Right'] = conns['index'].str.split(",").apply(lambda x: x[1])
        conns = conns.loc[:, ['Left', 'Right', 'Connection']]
        filtered = filter_significant_connections(conns, singletons)
        return filtered

    def get_itol_sig_above_connections(self, seq_type, sample, out_file=None,
                                       color='#ff0000', p_val=0.001, label='Label'):
        sig_abu_conns = self.get_significant_connections(seq_type, sample)
        sig_abu_conns = sig_abu_conns[sig_abu_conns['Lambda'] < sig_abu_conns['Observed']]
        sig_abu_conns = sig_abu_conns[sig_abu_conns['p-val'] <= p_val]
        sig_abu_conns['Observed'] = sig_abu_conns['Observed'].apply(str)
        sig_abu_conns['Color'] = color
        sig_abu_conns['Label'] = label
        sig_table = sig_abu_conns['Left'] + "," + sig_abu_conns['Right'] + "," +\
                    sig_abu_conns['Observed'] + "," + sig_abu_conns['Color'] +\
                    "," + sig_abu_conns['Label']
        sig_table = "\n".join(list(sig_table)).strip()
        conn_template = connection_template.format(label)
        sig_table = conn_template + sig_table
        if out_file:
            with open(out_file, 'w') as f:
                f.write(sig_table)
        else:
            return sig_table

    def get_itol_sig_below_connections(self, seq_type, sample, out_file=None,
                                       color='#0000ff', p_val=0.001, label='Label'):
        sig_abu_conns = self.get_significant_connections(seq_type, sample)
        sig_abu_conns = sig_abu_conns[sig_abu_conns['Lambda'] > sig_abu_conns['Observed']]
        sig_abu_conns = sig_abu_conns[sig_abu_conns['p-val'] <= p_val]
        sig_abu_conns['Observed'] = sig_abu_conns['Observed'].apply(str)
        sig_abu_conns['Color'] = color
        sig_abu_conns['Label'] = label
        sig_table = sig_abu_conns['Left'] + "," + sig_abu_conns['Right'] + "," +\
                    sig_abu_conns['Observed'] + "," + sig_abu_conns['Color'] +\
                    "," + sig_abu_conns['Label']
        sig_table = "\n".join(list(sig_table)).strip()
        conn_template = connection_template.format(label)
        sig_table = conn_template + sig_table
        if out_file:
            with open(out_file, 'w') as f:
                f.write(sig_table)
        else:
            return sig_table

    def get_popup_for_itol(self, out_file=None):
        tax = self.tax
        tax_table = tax['Kingdom'] + ";" + tax['Phylum'] + ";" + tax['Order'] + ";" +\
                    tax['Family'] + ";" + tax['Genus'] + ";" + tax['Species']
        tax_name = tax['Genus'].str.split(":").apply(lambda x: x[1]) + " " +\
                   tax['Species'].str.split(":").apply(lambda x: x[1].split(";")[0])
        popup_table = tax['OTU'] + "," + tax_name + "," + tax_table
        popup_str = "\n".join(list(popup_table)).strip()
        popup_str = popup_template + popup_str
        if out_file:
            with open(out_file, 'w') as f:
                f.write(popup_str)
        else:
            return popup_str


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
    with open(db_dir + 'HOMD_16S_rRNA_RefSeq_V14.5.qiime_spike.taxonomy', 'r') as t:
        for line in t:
            line = line.strip().split('\t')
            taxID = line[0]
            tax = line[1].strip().replace('__',':')
            tax = tax.replace(';',',')
            taxDict[taxID] = tax
    with open(db_dir + 'HOMD_16S_rRNA_RefSeq_V14.5.p9_spike.fasta', 'r') as f:
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


def process_unoise_fasta(input_file, seq_type):
    fst = io.read_fasta(input_file)
    acc = []
    for seq_id, seq in fst:
        seq_id_split = seq_id.split()
        bc_field = seq_id_split[-1]
        bc, otu = bc_field.split(";")[:2]
        bc = "barcode=" + bc.split("=")[1]
        new_seq_id = " ".join(seq_id_split[:-1]) + " " + bc + " sequence_type=" + seq_type + " OTU=" + otu
        acc.append([new_seq_id, seq])
    return acc
