#!/usr/bin/env python2
# coding: utf-8


import itertools
import re
import os, ssl
import urllib2, urllib
from collections import defaultdict
from bs4 import BeautifulSoup

if (not os.environ.get('PYTHONHTTPSVERIFY', '') and
    getattr(ssl, '_create_unverified_context', None)): 
    ssl._create_default_https_context = ssl._create_unverified_context




def GetCatalyticActivity(uniprots=None, infile=None):
	
	if uniprots:
		uniprots = uniprots.split(',')
		print uniprots
	elif infile:
		uniprots = loadFileasList(infile)
		print ' - Loading proteins from:' + infile


	accession2name = {}
	accession2descriptions = defaultdict(list)

	URL = "https://www.uniprot.org/uniprot/{}.xml"
	for uniprot in uniprots:	
		uniprot = uniprot.upper()
		url = URL.format( uniprot )
		print ' - Reading ' +  url	
		req = urllib2.Request(url)
		data = urllib2.urlopen(req).read()
		soup = BeautifulSoup(data, 'lxml')
		for entry in soup.find_all('entry'):
			if entry.gene is not None:
				# Extract descriptions including Bgee
				features = [f for f in entry.find_all() if f.has_attr('type') ]
				CA_feature = [f for f in features if f['type'] in ['catalytic activity'] ]
				if len(CA_feature) > 0:
					for ft in CA_feature:
						CAinfo = parseComment(ft)
						# Save accession
						accession =  str(entry.accession.get_text())
						# Save activeSite
						accession2descriptions[ accession ].append(CAinfo)
						# Save gene names
						names = [ str(g.get_text()) for g in entry.gene.find_all("name")]
						accession2name[ accession ] = names[0]
	return (accession2name,accession2descriptions)



def parseComment(feature):
    allReactionsInf = []

    if feature.find_all('reaction') is not None:
        for react in feature.find_all('reaction'):

            reactionInf = []
            # Extract reaction formula
            text = 'no reaction'
            if react.find('text') is not None:
                text = react.find('text').get_text()
            
            # Extract reaction rhea ID
            rheas = [f['id'] for f in react.find_all('dbreference') if f['type'] in 'Rhea']
            rhea_id = ';'.join(rheas)

            # Extract CHEBI ids of the participants
            participants = [f['id'] for f in react.find_all('dbreference') if f['type'] in 'ChEBI']
            part = ';'.join(participants)

            # Append reaction information
            reactionInf.append(text)
            reactionInf.append(rhea_id)
            reactionInf.append(part)
            allReactionsInf.append('\t'.join(reactionInf))
                
    return allReactionsInf

def loadFileasList(file):
	l = [line.strip() for line in open(file)]
	return l
