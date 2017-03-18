import unittest
import os
import pandas as pd
from .. import utilities, reads, io

os.chdir("epicBarcoder/tests/")

class TestUtilities(unittest.TestCase):
    def test_process_mapping(self):
        sampIDs, mapping, readCounts = utilities.process_mapping_file("test_map.txt")
        self.assertEqual(sampIDs, ['Sample1', 'Sample2', 'Sample3'])
        self.assertEqual(mapping, {'ACTTATTG': 'Sample1',
                                   'AATGTCAA': 'Sample2',
                                   'AGTATGCA': 'Sample3'})
        self.assertEqual(readCounts, {'ACTTATTG': 0,
                                      'AGTATGCA': 0,
                                      'AATGTCAA': 0})

    def test_process_fastq_and_mapping(self):
        utilities.process_fastq_and_mapping_file("test.fastq",
                                                 "output.fasta",
                                                 "test_map.txt",
                                                 "qual_summary.txt")
        with open("qual_summary.txt") as qo:
            qual_obs = list(qo)
        with open("summary.qual") as qc:
            qual_correct = list(qc)
        with open("output.fasta") as of:
            fasta_obs = list(of)
        with open("test.fasta") as tf:
            fasta_correct = list(tf)
        os.remove("qual_summary.txt")
        os.remove("output.fasta")
        self.assertEqual(qual_correct, qual_obs)
        self.assertEqual(fasta_correct, fasta_obs)

    def test_make_otus_and_assign(self):
        # Too time consuming for unit testing
        # Also assumes usearch in search path
        pass

    def test_process_unoise_fasta(self):
        processed_fasta_list = utilities.process_unoise_fasta("test_tax.fasta", '16S')
        first_entry = processed_fasta_list[0]
        header, seq = first_entry
        self.assertEqual(header.split()[-1], "OTU=Otu1")
        self.assertEqual(seq[:10], "CAGCAGCCGC")

    def test_fasta_to_bc_otu_table(self):
        processed_fasta_list = utilities.process_unoise_fasta("test_tax.fasta", '16S')
        io.write_fasta(processed_fasta_list, "processed.fasta")
        processed_table = utilities.fasta_to_bc_otu_table("processed.fasta")
        os.remove("processed.fasta")
        comp_list = [[">OM8s16_9", "OM8s16", "CCCCGGGGCTCCGGTGTCTG", "16S", "Otu1"],
                     [">OM8s16_10", "OM8s16", "AGGGGGGGCTCCGGTGTCTG", "16S", "Otu2"],
                     [">OM8s17_11", "OM8s17", "AAAATTTTCTCCGGTGTCTG", "16S", "Otu1"],
                     [">OM8s18_12", "OM8s18", "ATTCGAATGAAGGTTGCTTA", "16S", "Otu3"],
                     [">OM8s18_13", "OM8s18", "ATTCGAATGAAGGTTGCTTA", "16S", "Otu3"],
                     [">OM8s16_14", "OM8s16", "GGCCGGCCGGCCGTTGCTTA", "16S", "Otu1"],
                     [">OM8s16_15", "OM8s16", "GGCCGGCCGGCCGTTGCTTA", "16S", "Otu2"],
                     [">OM8s19_16", "OM8s19", "TTAATTAAATTAGTTGCTTA", "16S", "Otu4"],
                     [">OM8s19_17", "OM8s19", "TTAATTAAATTAGTTGCTTA", "16S", "Otu5"],
                     [">OM8s19_18", "OM8s19", "TTAATTAAATTAGTTGCTTA", "16S", "Otu6"],
                     [">OM8s19_19", "OM8s19", "TTAATTAAATTAGTTGCTTA", "16S", "Otu6"]]
        comp_tbl = pd.DataFrame(comp_list, columns=['Read', 'Sample', 'Barcode', 'Type', 'OTU'])
        self.assertEqual(comp_tbl.to_string(), processed_table.to_string())

    def test_get_grouped_table(self):
        processed_fasta_list = utilities.process_unoise_fasta("test_tax.fasta", '16S')
        io.write_fasta(processed_fasta_list, "processed.fasta")
        processed_table = utilities.fasta_to_bc_otu_table("processed.fasta")
        os.remove("processed.fasta")
        grouped_table = utilities.get_grouped_table(processed_table)
        items = list(grouped_table)
        self.assertEqual(items[0], ['Otu2'])
        self.assertEqual(set(items[2]), {'Otu1', 'Otu2'})
        self.assertEqual(items[3], ['Otu1'])
        self.assertEqual(set(items[4]), {'Otu3'})
        self.assertEqual(set(items[5]), {'Otu4', 'Otu5', 'Otu6'})

    def test_get_singletons(self):
        processed_fasta_list = utilities.process_unoise_fasta("test_tax.fasta", '16S')
        io.write_fasta(processed_fasta_list, "processed.fasta")
        processed_table = utilities.fasta_to_bc_otu_table("processed.fasta")
        os.remove("processed.fasta")
        grouped_table = utilities.get_grouped_table(processed_table)
        singletons = utilities.get_singletons(grouped_table)
        singleton_list = sorted(list(singletons['OTU']))
        self.assertEqual(singleton_list, ['Otu1', 'Otu1', 'Otu2'])

    def test_get_connections(self):
        processed_fasta_list = utilities.process_unoise_fasta("test_tax.fasta", '16S')
        io.write_fasta(processed_fasta_list, "processed.fasta")
        processed_table = utilities.fasta_to_bc_otu_table("processed.fasta")
        os.remove("processed.fasta")
        grouped_table = utilities.get_grouped_table(processed_table)
        connections = utilities.get_connections(grouped_table)
        self.assertEqual(list(connections['OTU']), [['Otu1', 'Otu2'],
                                                    ['Otu4', 'Otu5', 'Otu6']])


    def test_expand_connections(self):
        processed_fasta_list = utilities.process_unoise_fasta("test_tax.fasta", '16S')
        io.write_fasta(processed_fasta_list, "processed.fasta")
        processed_table = utilities.fasta_to_bc_otu_table("processed.fasta")
        os.remove("processed.fasta")
        grouped_table = utilities.get_grouped_table(processed_table)
        connections = utilities.get_connections(grouped_table)
        obs_conns = utilities.expand_connections(connections)
        exp_conns = [['OM8s16', 'Otu1,Otu2'],
                     ['OM8s19', 'Otu4,Otu5'],
                     ['OM8s19', 'Otu4,Otu6'],
                     ['OM8s19', 'Otu5,Otu6']]
        exp_conns = pd.DataFrame(exp_conns, columns=['Sample', 'Connection'])
        self.assertEqual(exp_conns.to_string(), obs_conns.to_string())

    def test_BarcodeContainer_unoise_singletons(self):
        container = utilities.BarcodeContainer(input_16S="test_tax.fasta", unoise=True)
        singleton_list = sorted(list(container.get_singletons('16S')['OTU']))
        self.assertEqual(singleton_list, ['Otu1', 'Otu1', 'Otu2'])

    def test_BarcodeContainer_unoise_connections(self):
        container = utilities.BarcodeContainer(input_16S="test_tax.fasta", unoise=True)
        connections = container.get_connections('16S')
        self.assertEqual(list(connections['Connection']), ['Otu1,Otu2',
                                                           'Otu4,Otu5',
                                                           'Otu4,Otu6',
                                                           'Otu5,Otu6'])

    def test_BarcodeContainer_get_samples(self):
        container = utilities.BarcodeContainer(input_16S="test_tax.fasta", unoise=True)
        samples = container.get_samples()
        self.assertEqual(samples, ['OM8s16', 'OM8s17', 'OM8s18', 'OM8s19'])

    def test_BarcodeContainer_get_itol_connections(self):
        container = utilities.BarcodeContainer(input_16S="test_tax.fasta", unoise=True)
        itol_connections = container.get_total_itol_connections('16S', 'OM8s19', '#ff0000')
        with open("test_itol_total_connections.txt") as f:
            test_itol_connections = ''.join(list(f)).strip()
        self.assertEqual(test_itol_connections, itol_connections)

    def test_BarcodeContainer_get_itol_abunds(self):
        container = utilities.BarcodeContainer(input_16S="test_tax.fasta", unoise=True)
        itol_abunds = container.get_itol_abunds('16S', 'OM8s16', '#ff0000')
        with open("test_itol_abund.txt") as f:
            test_itol_abunds = ''.join(list(f)).strip()
        self.assertEqual(test_itol_abunds, itol_abunds)

    def test_filter_significant_connections(self):
        conn = pd.read_csv("test_connections.csv", header=None)
        abu = pd.read_csv("test_singletons.csv", header=None)
        sig_conns = utilities.filter_significant_connections(conn, abu)
        self.assertAlmostEqual(sig_conns['p-val'][0], 5.9e-82)

    def test_BarcodeContainer_get_significant_connections(self):
        container = utilities.BarcodeContainer(input_16S="test_tax.fasta", unoise=True)
        sig_conns = container.get_significant_connections('16S', 'OM8s16')
        self.assertAlmostEqual(sig_conns['p-val'][0], 0.303265, places=5)


class TestReads(unittest.TestCase):
    def test_importFasta(self):
        imported_reads = reads.importFasta("test.fasta")
        first_read = imported_reads[0]
        self.assertEqual(first_read.header[:10], ">Sample1_0")
        self.assertEqual(first_read.seq[:10], "ACGCCACGGC")

    def test_filtBarcodePrimers(self):
        imported_reads = reads.importFasta("test.fasta")
        filt_reads = reads.filtBarcodePrimers(imported_reads, 20,
                                              'GATCATGACCCATTTGGAGAAGATG',
                                              'GGACTACHVGGGTWTCTAAT')
        self.assertEqual(filt_reads[0].header[:10], ">Sample2_2")
        self.assertEqual(filt_reads[0].seq[:10], "CAGCAGCCGC")


class TestIO(unittest.TestCase):
    def test_exportFasta(self):
        imported_reads = reads.importFasta("test.fasta")
        filt_reads = reads.filtBarcodePrimers(imported_reads, 20,
                                              'GATCATGACCCATTTGGAGAAGATG',
                                              'GGACTACHVGGGTWTCTAAT')
        io.exportFasta(filt_reads, "filtered.fasta")
        with open("filtered.fasta") as fo:
            fasta_filtered_obs = list(fo)
        with open("test_filtered.fasta") as tf:
            fasta_filtered_correct = list(tf)
        os.remove("filtered.fasta")
        self.assertEqual(fasta_filtered_obs, fasta_filtered_correct)
