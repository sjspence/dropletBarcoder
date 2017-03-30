from .reads import importFasta, removeFwdPrimer, removeFwdRevPrimer, \
	filtBarcodePrimers, trimLength, selectSamples, removeSamples, \
	splitByDegenerate
from .pairedEnds import pairConcatenate, readPairList
from .io import exportFasta, exportPairedFasta, exportOTUtable, \
	importQiimeOTU
from .utilities import clusterWithUsearch, filter_significant_connections, \
        process_fastq_and_mapping_file, add_otus_to_fasta, \
        output_abunds_and_connections, output_functions
#       move_barcodes_and_type_to_fasta_id
from .usearch_io import importClusterFast
from .taxonomy import importSintax, importTaxonomy
from .barcodes import createBarcodeDict, summarizeBarcoding
from .otuTables import buildOTUtable, invertHits
from .itol import itolHeatmap
from .trees import makeTreeConstraint, alignmentToSequence
from .dereplicate import getUniqueSeqs, expandDenoised, uniqueSeqsToOTU
from .parallel import run_array_job, run_batch_job
