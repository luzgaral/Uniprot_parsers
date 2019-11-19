#!/usr/bin/env python2
# coding: utf-8


import re




def PareseBgeeExpression(accession2name,accession2descriptions, outfileP):
	outf = open(outfileP, 'w')
	for prot in accession2descriptions:
		print ' - ' + prot
		for BgeeExpression in accession2descriptions[ prot ]:
			newlines = [prot] +  [accession2name[prot]] + [BgeeExpression]
			outf.write('\t'.join(newlines) + '\n')
	outf.closed

