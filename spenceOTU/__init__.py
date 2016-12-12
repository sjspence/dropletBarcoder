from .reads import importFasta, removeFwdPrimer, removeFwdRevPrimer, \
	trimLength, selectSamples, removeSamples, splitByDegenerate
from .pairedEnds import pairConcatenate, readPairList
from .io import exportFasta, exportPairedFasta, exportOTUtable
from .utilities import clusterWithUsearch
from .otuTables import importUsearchHits, importTaxonomy, buildOTUtable, \
	invertHits
from .itol import itolHeatmap
from .trees import makeTreeConstraint, alignmentToSequence
