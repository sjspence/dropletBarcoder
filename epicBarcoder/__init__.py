from .reads import importFasta, exportFasta, removeFwdPrimer, \
	removeFwdRevPrimer, filtBarcodePrimers, trimLength, selectSamples, \
	removeSamples, splitByDegenerate
from .pairedEnds import pairConcatenate, readPairList
from .io import exportPairedFasta, exportOTUtable, \
	importQiimeOTU, read_fasta, read_fastq, \
        write_fasta, write_fastq
from .utilities import clusterWithUsearch, filter_significant_connections, \
        process_fastq_and_mapping_file, add_otus_to_fasta, \
        output_functions, fasta_to_bc_otu_table, \
        get_grouped_table, make_otus_and_assign, \
        process_unoise_fasta, BarcodeContainer
from .otuTables import buildOTUtable, invertHits
from .usearch_io import importClusterFast
from .itol import itolHeatmap, itolSimpleBar, itolConnections, itolHover
from .trees import makeTreeConstraint, alignmentToSequence, tOTU_pickRepSeqs
from .dereplicate import getUniqueSeqs, uniqueSeqsToOTU, otuToHeaders, \
        expandDenoised
from .taxonomy import importSintax, tOTUmap, importTaxonomy
from .barcodes import createBarcodeDict, summarizeBarcoding, \
	tOTU_singletonAbundances, tOTU_quantifyPairs, pickSigPairs
from .parallel import run_array_job, run_batch_job
