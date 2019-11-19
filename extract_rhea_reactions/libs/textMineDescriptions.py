#!/usr/bin/env python2
# coding: utf-8


import re




def PareseRhea(accession2name,accession2descriptions, outfileP):
	outf = open(outfileP, 'w')
	for prot in accession2descriptions:
		print ' - ' + prot
		for reaction in accession2descriptions[ prot ]:
			newlines = [prot] +  [accession2name[prot]] + reaction
			outf.write('\t'.join(newlines) + '\n')
	outf.closed

