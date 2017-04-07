#!/usr/bin/python
from ete3 import NCBITaxa
import sys
from tqdm import tqdm
import gzip
ncbi=NCBITaxa()

def returnLineages(taxid):
	return ncbi.get_lineage(int(taxid))

def lastCommon(tax1,tax2):
	return [t1 for t1,t2 in zip(tax1,tax2) if t1==t2]

def lca(reads):
	# get list of taxes from centrifuge reads
	taxes=[r.split('\t')[2] for r in reads if r.split('\t')[2] != '0']
	if len(taxes)<=1:
		try:
			return taxes[0]
		except: return 0

	# get lineages from taxids
	lins=map(returnLineages,taxes)

	# get last common acestor of first two lineages
	pl=reduce(lastCommon, lins)
	return str(pl[-1])

def processFile(f):
	reads=gzip.open(f,'r').read()
	reads=reads.split('\n')[1:-1]
	preName = None
	for i in tqdm(xrange(len(reads))):
		c=reads[i].split('\t')
		if int(c[-1]) > 1 and c[0] != preName:
			preName=c[0]
			c[2]=lca(reads[i:i+int(c[-1])])
			c[-1]="1"
			print "\t".join(c)
		elif c[0] != preName:
			print reads[i]

processFile(sys.argv[1]) 
