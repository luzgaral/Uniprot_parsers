#!/usr/bin/env python2
# coding: utf-8


import itertools
import re
from collections import defaultdict
import urllib2, urllib
from bs4 import BeautifulSoup
from libs.textMineDescriptions import loadFileasList




def GetMutationsData(outFile, infile=None):
	uniprots = loadFileasList(infile)
	print ' - Loading proteins from:' + infile


	name2accession = {}
	accession2name = {}
	accession2descriptions = defaultdict(list)

	outf = open(outFile, 'w')
	header = ['#uniprot_accession', 'gene_symbol', 'variant_type' , 'description', 'PubMed_evidence', 'wt_allele', 'alt_allele', 'start_position', 'end_position']
	outf.write('\t'.join(header) + '\n')
	URL = "http://www.uniprot.org/uniprot/{}.xml"
	for uniprot in uniprots:	
		uniprot = uniprot.upper()
		url = URL.format( uniprot )
		print ' - Reading ' +  url	
		req = urllib2.Request(url)
		data = urllib2.urlopen(req).read()
		soup = BeautifulSoup(data, "html.parser")
		for entry in soup.find_all('entry'):
			if entry.gene is not None:
				# Get accession
				accession =  str(entry.accession.get_text())
				# Get gene names and synonyms
				names = [ str(g.get_text()) for g in entry.gene.find_all("name")]
				for n in names:
					name2accession[ n ] = accession
				accession2name[ accession ] = names[0]
				# Build PubMed evidences dictionary
				evidenceDic = parseEvidences(entry)
				# Extract mutation and variant features
				features = [f for f in entry.find_all('feature') if f['type'] in ['sequence variant', 'mutagenesis site'] ] 
				for ft in features:
					description = parseMutFeature(ft, evidenceDic)
					if description is not None:
						accession2descriptions[ accession ].append(description)
						out = [accession, names[0], description['type'] , description['description'], description['evidence'], description['wt'], ';'.join(description['alt']), description['start'], description['end']]
						outf.write('\t'.join(out) + '\n')
	outf.closed
	return (name2accession,accession2name,accession2descriptions)



def parseMutFeature(feature, evidenceDic):
	
	description = None
	
	if feature.has_attr('description') and feature.find('original') is not None:

		desc = feature['description']
		wt = feature.find('original').get_text()
		alt = [ str(alt.get_text()) for alt in feature.findAll('variation')]
		
		evidence = None
		if feature.has_attr('evidence'):
			evidence = [ evidenceDic[ k ] for k in feature['evidence'].split() if evidenceDic[ k ]]
			evidence = list(itertools.chain(*evidence))
			evidence = ';'.join(evidence)
		
		if feature.find('position') is not None:
			start = feature.find('position')['position']
			end = start
		elif feature.find('begin') is not None and feature.find('end') is not None:
			start = feature.find('begin')['position']
			end = feature.find('end')['position']
		
		description = {'type':str(feature['type']), 
					   'description': str(desc),
					   'evidence':str(evidence),
					   'wt':str(wt), 
					   'alt':alt, 
					   'start':str(start), 'end':str(end)}
	return description


def parseEvidences(entry):
	evidenceDic = {}
	for e in entry.find_all('evidence'):
		ekey = e['key']
		pubmedList = []
		for s in e.find_all('source'):
			pubmedList =  [ str(p['id']) for p in s.find_all('dbreference') if p['type'] in 'PubMed']
		evidenceDic[ ekey ] = pubmedList
	return evidenceDic
