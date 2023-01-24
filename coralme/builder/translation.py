import Bio
import coralme

#from functools import reduce
#def reducer(accumulator, element):
	#for key, value in element.items():
		#accumulator[key] = accumulator.get(key, 0) + value
	#return accumulator

from collections import Counter

def add_subreactions_to_model(me_model, subreactions):
	# add subreactions associated with translation initiation, elongation, termination and postprocessing
	for subreaction in subreactions:
		for rxn, info in subreaction.items():
			data = coralme.core.processdata.SubreactionData(rxn, me_model)
			data.enzyme = info['enzymes']
			if me_model.global_info['translation_subreactions'][rxn] == '':
				data.stoichiometry = {}
			else:
				data.stoichiometry = me_model.process_data.get_by_id(me_model.global_info['translation_subreactions'][rxn]).stoichiometry
				data._element_contribution = data.calculate_element_contribution() #info.get('element_contribution', {})

	#old code
	#for rxn, info in termination_subreactions.items():
		#data = coralme.core.processdata.SubreactionData(rxn, me_model)
		#data.enzyme = info['enzymes']
		#data.stoichiometry = info['stoich']
		#data._element_contribution = info.get('element_contribution', {})

def add_charged_trna_subreactions(me_model, organelle = 'c', transl_table = set([11]), translation_stop_dict = {}, special_trna_subreactions = {}):
	"""
	Create subreaction for each codon. this will be used to model
	the addition of charged tRNAs to the elongating peptide
	"""
	#for codon in coralme.util.dogma.codon_table:
	for codon in me_model.global_info['stop_codons']:
		#if coralme.util.dogma.codon_table[codon] == '*':
		stop_codon = codon.replace('T', 'U')
		stop_enzyme = translation_stop_dict.get(stop_codon, 'CPLX_dummy')
		me_model.add_metabolites([coralme.core.component.Complex(stop_enzyme)])
		subreaction_data = coralme.core.processdata.SubreactionData(stop_codon + '_' + stop_enzyme + '_mediated_termination_' + organelle, me_model)
		subreaction_data.enzyme = stop_enzyme
		subreaction_data.stoichiometry = {}
		#else:

	codon_table = Bio.Data.CodonTable.generic_by_id[list(transl_table)[0]]

	#for codon, aa in {k:v for k,v in me_model.global_info['codon_table'].forward_table.items() if 'U' not in k}.items():
	for codon, aa in {k:v for k,v in codon_table.forward_table.items() if 'U' not in k}.items():
		#full_aa = coralme.util.dogma.amino_acids[coralme.util.dogma.codon_table[codon]]
		full_aa = coralme.util.dogma.amino_acids[aa]

		if not me_model.metabolites.has_id(full_aa + '_' + organelle):
			continue

		#aa = full_aa.split('_')[0]
		subreaction_data = coralme.core.processdata.SubreactionData(full_aa + '_' + organelle + '_addition_at_' + codon.replace('T', 'U'), me_model)
		trna = 'generic_tRNA_' + codon.replace('T', 'U') + '_' + full_aa + '_' + organelle
		# Default AA loader enzyme
		subreaction_data.enzyme = me_model.global_info['amino_acid_loader']
		# Accounts for GTP hydrolyzed by EF-TU and the ATP hydrolysis to AMP required to add the amino acid to the tRNA
		#subreaction_data.stoichiometry = reduce(
			#reducer, [{trna : -1}, me_model.global_info['atp_trna_loading'], me_model.global_info['gtp_hydrolysis']], {})
		subreaction_data.stoichiometry = Counter({trna : -1})
		subreaction_data.stoichiometry.update(me_model.process_data.get_by_id('atp_hydrolysis_trna_loading').stoichiometry)
		subreaction_data.stoichiometry.update(me_model.process_data.get_by_id('gtp_hydrolysis').stoichiometry)
		subreaction_data._element_contribution = subreaction_data.calculate_element_contribution()

	# Add subreactions for start codon and selenocysteine
	for rxn, info in special_trna_subreactions.items():
		data = coralme.core.processdata.SubreactionData(rxn, me_model)
		data.enzyme = info['enzymes']
		data.stoichiometry = info['stoich']
		data._element_contribution = info.get('element_contribution', {})
