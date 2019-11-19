#!/usr/bin/env python2
# coding: utf-8


import re




def ParsePTMsDescriptions(accession2name,accession2descriptions, outfileP):
	outf = open(outfileP, 'w')
	for prot in accession2descriptions:
		print ' - ' + prot
		newlines = []
		for text in accession2descriptions[ prot ]:
			# ----- Get descriptions to text mining
			descriptionList = get_descriptionList(text)
			# ----- Analyse each descritpion separately
			for d in descriptionList:
				data = parseData(d)
				for ptm in data:
					for site in data[ptm]:
						n = [prot] + [accession2name[prot]] + [site]+ [PTMs[ptm]]
						n = '\t'.join(n)
						newlines.append(n)
		for l in list(set(newlines)):
			outf.write(l + '\n')
	outf.closed




def parseData(desc):
	data = {}
	for ptm in PTMs:
		if re.search(ptm, desc, re.IGNORECASE):
			desc = re.sub('/', ' ', desc)
			words = desc.split(' ')
			sites = [s for s in words if re.search(r'[A-Z][a-z][a-z]-[1-9]', s)]
			sites = [re.sub(r',', '', s) for s in sites]
			sites = [re.sub(r';', '', s) for s in sites]
			sites = [re.sub(r'\.', '', s) for s in sites]
			sites = [re.sub(r'\"', '', s) for s in sites]
			sites = [re.sub(r'\'', '', s) for s in sites]
			sites = [re.sub(r'\)', '', s) for s in sites]
			sites = [re.sub(r'\(', '', s) for s in sites]
			sites = ['-'.join(s.split('-')[0:2]) for s in sites]
			for a3 in AA:
				sites = [re.sub(a3+'-', AA[a3]+'\t', s) for s in sites]
			sites = [s for s in sites if not re.search(r'-', s)]
			if sites is not None:
				data[ptm] = list(set(sites))
	return data




def get_descriptionList(st):
	st = re.sub('promote', '. ', st)
	descriptionList = st.split('. ')
	descriptionList = [s for s in descriptionList if re.search(r'-[1-9]', s)]
	return descriptionList



def loadFileasList(file):
	l = [line.strip() for line in open(file)]
	return l


def loadFileasDic(file):
	dic = {}
	with open(file) as f:
		for line in f:
			s = line.strip().split('\t')
			dic[s[0]] = s[1]
	return dic

PTMs = loadFileasDic('data/ptm.annot')
AA = loadFileasDic('data/aa.annot')

