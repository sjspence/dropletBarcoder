from .reads import importFasta, removeFwdPrimer, removeFwdRevPrimer, \
	trimLength, selectReads
from .pairedEnds import pairConcatenate, readPairList
from .io import exportFasta, exportPairedFasta
from .utilities import clusterWithUsearch
from .otuTables import importUsearchHits, importTaxonomy, buildOTUtable, \
	invertHits
from .itol import itolHeatmap
from .trees import makeTreeConstraint
