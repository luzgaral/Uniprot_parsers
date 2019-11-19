#!/usr/bin/env python2
# coding: utf-8


import re
from collections import OrderedDict





def ParseMutInfo(name2accession,accession2name,accession2descriptions, outfileP):

	outf = open(outfileP, 'w')
	header = ['#uniprot_accession', 'gene_symbol', 'variant_type' , 'start_position', 'end_position', 'wt_allele', 'alt_allele', 'regex', 'effect_sign', 'target_uniprot_accession', 'target_name', 'cleaned_description', 'PubMed_evidence']
	outf.write('\t'.join(header) + '\n')

	for prot in accession2descriptions:
		
		print ' - ' + prot

		for data in accession2descriptions[ prot ]:
			
			# ----- Protein and Description info
			pubmed = data['evidence']
			new_aa = data['alt']
			info = [prot, accession2name[prot], data['type'], data['start'], data['end'], data['wt']]

			# ----- Get descriptions to text mining
			descriptionList = get_descriptionList(data['description'], accession2name[prot])

			# ----- Analyse each descritpion separately
			for d in descriptionList:

				## ----- Unify synonyms
				d = unifySynonyms(d, accession2name[prot])
				
				## ----- STABILITY
				if 'stabil' in d or 'stable' in d or 'protein processing' in d:
					# KNOWN side EFFECT (decrease)
					for reg in regex_stab:
						if re.search(reg, d, re.IGNORECASE) and not re.search('not .*'+reg, d, re.IGNORECASE) and not re.search('no .*'+reg, d, re.IGNORECASE) and not re.search('compensates for '+reg, d, re.IGNORECASE):
							for new in new_aa:
								newline = info + [new] + [reg] + ['destabilizing'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')
					for reg in regex_stab_unknown:
						# UNKNOWN side EFFECT
						if re.search(reg, d, re.IGNORECASE) and not re.search('not '+reg, d, re.IGNORECASE) and not re.search('no '+reg, d, re.IGNORECASE):
							for new in new_aa:
								newline = info + [new] + [reg] + ['stability_effect'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')
						## NO EFFECT
						if re.search('not '+reg, d, re.IGNORECASE) or re.search('no '+reg, d, re.IGNORECASE) or re.search('no apparent '+reg, d, re.IGNORECASE):
							for new in new_aa:
								newline = info + [new] + [reg] + ['no_effect'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')


				## ----- ACTIVITY
				if ('activ' in d or ' function' in d or 'suppressive ability' in d or 'mutant protein is more efficient' in d) and not 'retain' in d:
					#print d
					# LOSS
					for reg in regex_loss_actv:
						if (' activ' in d or 'function' in d or 'suppressive ability' in d) and re.search(reg, d, re.IGNORECASE) and not re.search('not .*'+reg, d, re.IGNORECASE) and not re.search('no .*'+reg, d, re.IGNORECASE) :
							for new in new_aa:
								newline = info + [new] + [reg + ' activity' ] + ['loss'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')
					# GAIN
					for reg in regex_gain_actv:
						if re.search(reg, d, re.IGNORECASE) and not any(re.search(ele, d, re.IGNORECASE) for ele in confGain) and not re.search('not .*'+reg, d, re.IGNORECASE) and not re.search('no .*'+reg, d, re.IGNORECASE):
							for new in new_aa:
								newline = info + [new] + [reg ] + ['gain'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')
					# NO EFFECT
					for reg in regex_loss_actv + regex_gain_actv :
						if re.search('not '+reg, d, re.IGNORECASE) or re.search('no .*'+reg, d, re.IGNORECASE):
							for new in new_aa:
								newline = info + [new] + [reg + ' activity'] + ['no_effect'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')
					if re.search(r'alter([\sA-Za-z]*) activity', d, re.IGNORECASE) and not re.search(r'no alter([\sA-Za-z]*) activity', d, re.IGNORECASE) and not re.search(r'not alter([\sA-Za-z\-]*) activity', d, re.IGNORECASE):
						for new in new_aa:
							newline = info + [new] + ['alters activity'] + ['effect'] +['null'] + ['null']  + [d] + [pubmed]
							outf.write('\t'.join(newline) + '\n')
				if 'inactiv' in d:
					# GAIN inactivity
					for reg in regex_gain_inactv:
						if re.search(reg, d, re.IGNORECASE):
							for new in new_aa:
								newline = info + [new] + [reg] + ['loss'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')
					# LOSS inactivation
					for reg in regex_loss_inactv:
						if re.search(reg, d, re.IGNORECASE):
							for new in new_aa:
								newline = info + [new] + [reg] + ['gain'] +['null'] + ['null']  + [d] + [pubmed]
								outf.write('\t'.join(newline) + '\n')



				## ----- BINDING and PHOSPO
				if ('bind' in d or ' interact' in d or ' phosphorylat' in d) and ' between ' not in d and not 'phosphorylation of ' in d:
					# search for any regex expression
					# LOSS
					for reg in regexp_loss:
						if re.search(reg, d, re.IGNORECASE) and not re.search('not .*'+reg, d, re.IGNORECASE) and not re.search('no .*'+reg, d, re.IGNORECASE):
							newlines = parseData_with_interactors(reg, d, new_aa, info, pubmed, name2accession, effect='loss')
							if len(newlines) > 0:	
								outf.write('\n'.join(newlines) + '\n')



					# GAIN
					for reg in regexp_gain:
						if not any(re.search(ele, d, re.IGNORECASE) for ele in confGain) and re.search(reg, d, re.IGNORECASE) and not re.search('not .*'+reg, d, re.IGNORECASE) and not re.search('no .*'+reg, d, re.IGNORECASE):
							newlines = parseData_with_interactors(reg, d, new_aa, info, pubmed, name2accession, effect='gain')
							if len(newlines) > 0:	
								outf.write('\n'.join(newlines) + '\n')

					# NO EFFECT
					for reg in regexp_gain + regexp_loss:
						if re.search('not .*'+reg, d, re.IGNORECASE) or re.search('no .*'+reg, d, re.IGNORECASE):
							newlines = parseData_with_interactors(reg, d, new_aa, info, pubmed, name2accession, effect='no_effect')
							if len(newlines) > 0:	
								outf.write('\n'.join(newlines) + '\n')

	outf.closed




def remove_lower(lst):
	return [w for w in lst if w.isupper()]






def getInteractors(interactors):
	interactors = re.sub(r'histone ', '', interactors)
	interactors = re.sub(r'-[0-9] ', ' ', interactors)
	interactors = interactors.split('-binding')[0]
	interactors = interactors.split('-mediated')[0]
	interactors = interactors.split('-induced')[0]
	interactors = interactors.split('lead')[0]
	interactors = interactors.split('not')[0]
	interactors = interactors.split('to the mutant')[0]
	interactors = interactors.split('in response to')[0]
	interactors = re.findall(r"[\w']+", interactors)
	interactors = [s for s in interactors if len(s) >= 3]
	interactors = list(OrderedDict.fromkeys(remove_lower(interactors)))
	return list(set(interactors))






def parseData_with_interactors(reg, d, new_aa, info, pubmed, name2accession, effect='no_effect', add=''):
	# identify interactor poteins or enzymes in the description
	if len(reg.split(' ')) > 1 and re.search(r'[A-Z,0-9] '+reg.split(' ')[1], d):
		string = re.compile(reg.split(' ')[1], re.IGNORECASE).split(d)
		interactors = getInteractors(string[0])
	else:
		string = re.compile(reg, re.IGNORECASE).split(d)
		interactors = getInteractors(string[1])
	interactors = list(OrderedDict.fromkeys(interactors))
	interactors = filter(lambda i: not re.search(i, info[1]) , interactors)
	newlines = []
	# if there are interactors, extract and annotate
	for i in interactors:
		uniprot = 'unknown'
		if i in name2accession:
			uniprot = name2accession[i]
		for new in new_aa:
			newline = info + [new] + [reg + add] + [effect] +[uniprot] + [i]  + [d] + [pubmed]
			newlines.append('\t'.join(newline))
	# if no interactors, print interactors equal to null
	if len(interactors) == 0:
		for new in new_aa:
			newline = info + [new] + [reg + add] + [effect] + ['null'] + ['null']  + [d] + [pubmed]
			newlines.append('\t'.join(newline))
	return list(set(newlines))





def unifySynonyms(st, n):
	# ----- Unify synonyms
	st = re.sub('complex with', 'interaction with', st)
	st = re.sub('association to', 'interaction with', st)
	st = re.sub('association with', 'interaction with', st)
	st = re.sub('heterodimerization with', 'interaction with', st)
	st = re.sub('dimerization with', 'interaction with', st)
	st = re.sub('assembly with', 'interaction with', st)
   	st = re.sub('ability to bind', 'binding to', st)
	st = re.sub('affinity to', 'binding to', st)
	st = re.sub('binding affinity for', 'binding to', st)
	st = re.sub('interaction between ' + n + ' and ', 'interaction with ', st)
	st = re.sub('-mediated phosphorylation ', ' phosphorylation ', st)
	st = re.sub('-induced phosphorylation ', ' phosphorylation ', st)
	st = re.sub('-binding ', ' binding ', st)
	st = re.sub(' change in ', ' effect in ', st)
	st = re.sub('affinity for', 'binding to', st)
	st = re.sub(r'activation of ([A-Z])', r'interaction with \1', st)
	st = re.sub('binding activity', 'binding to', st)
	st = re.sub('binding activity', 'binding to', st)
	st = re.sub('regulation by', 'interaction with', st)
	st = re.sub('mildly', 'marginal', st)
	st = re.sub('slightly', 'marginal', st)
	st = re.sub('without', 'no', st)
	st = re.sub(r' nor ', r' no ', st)
	# Special cases
	st = re.sub('lipid and protein', 'lipid-and-protein', st)
	return st








def get_descriptionList(st, prot): # DIRTY - Debug mode
	## ----- Unify synonyms
	st = unifySynonyms(st, prot)
	st = st.replace(', while it', '. ')
	# ----- split text into independent descriptions 
	# ----- independent descriptions are separated mostly by "."
	# ----- but also by
	#	   "and" and ", " (if not surrounded by genes)	
	st = re.sub(r', ([a-z])', r'. \1', st)
	st = re.sub(r' and ([a-z])', r'. \1', st)
	st = re.sub(r' but ([a-z])', r'. \1', st)
	#	   ";" (no by "; when")
	st = st.replace('; when', ' when').replace('; ', '. ').replace('leading to ', '. ')
	descriptionList = st.split('. ')
	return descriptionList



def loadFileasList(file):
	l = [line.strip() for line in open(file)]
	return l







## LOSS edge
regexp_loss = loadFileasList('data/regex/ppi_loss')
## GAIN edge
regexp_gain = loadFileasList('data/regex/ppi_gain')
## OVERALL STABILITY
regex_stab = loadFileasList('data/regex/stab')
## UNKNOWN EFFECT OVER STABILITY
regex_stab_unknown = loadFileasList('data/regex/stab_unknown')
## LOSS activity
regex_loss_actv = loadFileasList('data/regex/activ_loss')
regex_gain_inactv = loadFileasList('data/regex/incative')
## Gain activity
regex_gain_actv = loadFileasList('data/regex/activ_gain')
regex_loss_inactv = loadFileasList('data/regex/incative_loss')
## General loss
general_loss = loadFileasList('data/regex/general_loss')
## Confounding
confounding = loadFileasList('data/regex/confounding')
confGain = general_loss + confounding
## General verbs
effect_verbs = loadFileasList('data/regex/effect_verbs')

