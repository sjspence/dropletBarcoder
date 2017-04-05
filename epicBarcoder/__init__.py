from .reads import importFasta, removeFwdPrimer, removeFwdRevPrimer, \
	filtBarcodePrimers, trimLength, selectSamples, removeSamples, \
	splitByDegenerate
from .pairedEnds import pairConcatenate, readPairList
from .io import exportFasta, exportPairedFasta, exportOTUtable, \
	importSintax, importQiimeOTU, read_fasta, read_fastq, \
        write_fasta, write_fastq
from .utilities import clusterWithUsearch, filter_significant_connections, \
        process_fastq_and_mapping_file, add_otus_to_fasta, \
        output_functions, fasta_to_bc_otu_table, \
        get_grouped_table, make_otus_and_assign, \
        process_unoise_fasta, BarcodeContainer
from .otuTables import importUsearchHits, importTaxonomy, buildOTUtable, \
	invertHits
from .itol import itolHeatmap
from .trees import makeTreeConstraint, alignmentToSequence
from .dereplicate import getUniqueSeqs, uniqueSeqsToOTU, otuToHeaders
from .parallel import run_array_job, run_batch_job
