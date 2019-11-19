#!/usr/bin/env python2
# coding: utf-8


import itertools
import re
import urllib2, urllib
from collections import defaultdict
from bs4 import BeautifulSoup
from libs.textMineDescriptions import loadFileasList,loadFileasDic




def GetPTMsDsescriptions(uniprots=None, infile=None):
	
	if uniprots:
		uniprots = uniprots.split(',')
		print uniprots
	elif infile:
		uniprots = loadFileasList(infile)
		print ' - Loading proteins from:' + infile


	accession2name = {}
	accession2descriptions = defaultdict(list)

	URL = "http://www.uniprot.org/uniprot/{}.xml"
	for uniprot in uniprots:	
		uniprot = uniprot.upper()
		url = URL.format( uniprot )
		print ' - Reading ' +  url	
		req = urllib2.Request(url)
		data = urllib2.urlopen(req).read()
		soup = BeautifulSoup(data)
		for entry in soup.find_all('entry'):
			if entry.gene is not None:
				# Extract descriptions incuding PTMs
				features = [f for f in entry.find_all('comment') if f['type'] in ['enzyme regulation', 'PTM'] ]
				if len(features) > 0:
					print '   Number of features: ' + str(len(features))
					for ft in features:
						description = parseComment(ft)
						if description is not None:
							# Save accession
							accession =  str(entry.accession.get_text())
							# Save description
							accession2descriptions[ accession ].append(description)
							# Save gene names
							names = [ str(g.get_text()) for g in entry.gene.find_all("name")]
							accession2name[ accession ] = names[0]
	return (accession2name,accession2descriptions)


def parseComment(feature):
	description = None
	if feature.find('text') is not None:
		text = feature.find('text').get_text()
		if re.search(r'[A-Z][a-z][a-z]-[1-9]', text) and any(re.search(ptm, text, re.IGNORECASE) for ptm in PTMs):
			description = text
	return description


PTMs = loadFileasDic('data/ptm.annot')