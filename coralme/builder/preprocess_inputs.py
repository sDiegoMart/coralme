#!/usr/bin/python3
import numpy
import pathlib
import pandas
import warnings
import logging
from Bio import SeqIO

try:
	warnings.simplefilter(action = 'ignore', category = pandas.errors.SettingWithCopyWarning)
except:
	warnings.warn("This pandas version does not allow for correct warning handling. Pandas >=1.5.1 is suggested.")

def generate_organism_specific_matrix(genbank, locus_tag, model):
	contigs = []
	for contig in SeqIO.parse(genbank, 'genbank'):
		contigs.append(contig)

	# get all features
	lst = [ x for y in [ x.features for x in contigs ] for x in y ]
	lst = [ x for x in lst if x.type in [ 'CDS', 'rRNA', 'tRNA', 'ncRNA', 'tmRNA', 'misc_RNA' ] ]

	# create a pandas DataFrame with organism-specific information to be completed with the builder data
	df = pandas.DataFrame(columns = [
		'Gene Locus ID',
		'Gene Names',
		'Old Locus Tag',
		'BioCyc',
		'Reference BBH',
		'Definition',
		'Feature Type',
		'Complex Name',
		'Complex ID',
		'Cofactors in Modified Complex',
		'Generic Complex ID',
		'MetaComplex ID',
		'ME-model SubReaction',
		'M-model Reaction ID',
		'Reaction Name',
		'Reversibility',
		'GroEL_dependent_folding',
		'DnaK_dependent_folding',
		'N_terminal_methionine_cleavage',
		'RNA mods/enzyme',
		'Complex Location',
		'Subunit Location',
		'Translocation Pathway',
		#'Translocation Multiplier'
	])

	def get_reaction(x):
		if x is None:
			return []
		else:
			lst = []
			for gene in x.split(';'):
				if model.genes.has_id(gene):
					lst.append([ x.id for x in model.genes.get_by_id(gene).reactions ])

			# (An)A (old)locustag/gene name can be associated to many reactions
			lst = [ x for y in lst for x in y ]
			if lst == []:
				return []
			else:
				return lst

	def get_reaction_name(x):
		if x is not None:
			try:
				return model.reactions.get_by_id(x).name
			except:
				return None
		else:
			return None

	def get_reversibility(x):
		if x is not None:
			try:
				return 'TRUE' if model.reactions.get_by_id(x).reversibility else None
			except:
				return None
		else:
			return None

	def get_spontaneous(x):
		if x is not None:
			if 'spontaneous' in model.reactions.get_by_id(x).name:
				return 'TRUE'
			else:
				return None
		else:
			return None

	#df['Gene Locus ID'] = [ x.qualifiers.get('locus_tag', [None])[0] for x in lst ]
	df['Gene Locus ID'] = [ x.qualifiers.get(locus_tag, [None])[0] for x in lst ]
	df['Definition'] = [ x.qualifiers.get('product', [None])[0] for x in lst ]
	df['Feature Type'] = [ x.type if x.qualifiers.get('pseudo') is None else 'pseudo' for x in lst ]

	tmp = [ x.qualifiers.get('gene', None) for x in lst ]
	df['Gene Names'] = [ ';'.join(x) if x is not None else None for x in tmp ]

	tmp = [ x.qualifiers['old_locus_tag'] if x.qualifiers.get('old_locus_tag', None) is not None else None for x in lst ]
	df['Old Locus Tag'] = [ ';'.join(x) if x is not None else None for x in tmp ]

	df['M-model Reaction ID'] = df['Old Locus Tag'].apply(lambda x: get_reaction(x))
	df['M-model Reaction ID'] += df['Gene Names'].apply(lambda x: get_reaction(x))
	df['M-model Reaction ID'] += df['Gene Locus ID'].apply(lambda x: get_reaction(x))
	df = df.explode('M-model Reaction ID')

	df['Reaction Name'] = df['M-model Reaction ID'].apply(lambda x: get_reaction_name(x))
	df['Reversibility'] = df['M-model Reaction ID'].apply(lambda x: get_reversibility(x))

	# df.set_index(['Gene Locus ID', 'Definition', 'Feature type'], inplace = True)
	return df.sort_values(['M-model Reaction ID', 'Gene Locus ID'])

def complete_organism_specific_matrix(builder, data, model, output):
	# ME-model homology to reference
	def bbh(x, dct, keys):
		#tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = x[keys].to_list()
		tags = [ str(x).split(';') for x in tags ]

		lst = []
		for tag in [ x for y in tags for x in y ]:
			lst.append(dct.get(tag, None))
		lst = [ x for x in lst if x is not None ]

		if len(lst) != 0:
			return lst[0]

	if builder.configuration.get('biocyc.genes', False):
		# We reuse the bbh function, but changed the dictionary of relationships between IDs
		dct = { v['Accession-1']:k for k,v in builder.org.gene_dictionary.iterrows() }
		data['BioCyc'] = data.apply(lambda x: bbh(x, dct, keys = ['Gene Locus ID', 'Old Locus Tag', 'Gene Names']), axis = 1)
	else:
		data['BioCyc'] = None

	if hasattr(builder, 'homology'):
		dct = builder.homology.mutual_hits
		data['Reference BBH'] = data.apply(lambda x: bbh(x, dct, keys = ['Gene Locus ID', 'Old Locus Tag', 'BioCyc']), axis = 1)

	# ME-model complexes: restructure complexes_df to obtain the correct complex stoichiometry from the index
	df = builder.org.complexes_df.copy(deep = True)
	#df = df[~df.index.str.contains('MONOMER')]
	df = df[df['genes'].str.contains('\(')]
	df['genes'] = df['genes'].str.split(' AND ')
	df = df.explode('genes')
	df['stoich'] = df['genes'].apply(lambda x: '1' if x.split('(')[1][:-1] == '' else str(x.split('(')[1][:-1]))
	df['genes'] = df['genes'].apply(lambda x: x.split('(')[0].replace('-MONOMER', ''))
	df.index = df.index + ':' + df['stoich']

	def complexes(x, df):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = [ str(x).split(';') for x in tags ]

		cplx = []
		for tag in [ x for y in tags for x in y ]:
			res = df[df['genes'].str.fullmatch(tag)].index.to_list()
			if len(res) != 0:
				cplx.append(res)
		if len(cplx) != 0:
			return cplx[0]

	data['Complex ID'] = data.apply(lambda x: complexes(x, df), axis = 1)
	data = data.explode('Complex ID')

	# ME-model cofactors
	def cofactors(x, dct):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], str(x['Complex ID']).split(':')[0] ]
		tags = [ str(x).split(';') for x in tags ]

		cplxs = []
		for tag in [ x for y in tags for x in y ]:
			cplxs.append(dct.get(tag, None))

		cplxs = [ x for x in cplxs if x is not None ]
		if len(cplxs) != 0:
			mods = [ x.split('_mod_')[1:] for x in cplxs ][0]
			mods = [ x for x in mods if not x.startswith('3hocta') ] # metabolic modification in ACP
			mods = [ x for x in mods if not x.startswith('Oxidized') ] # metabolic modification in ferredoxin and other proteins
			mods = [ x for x in mods if not x.startswith('palmitate') ] # metabolic modification from 2agpg160 in the lpp gene
			mods = [ x.replace('lipo', 'lipoyl') for x in mods ]
			logging.warning('The modification \'lipo\' was renamed to \'lipoyl\'.')
			mods = [ x.replace('NiFeCoCN2', 'NiFe_cofactor') for x in mods ]
			logging.warning('The modification \'NiFeCoCN2\' was renamed to \'NiFe_cofactor\'.')
			mods = [ '{:s}(1)'.format(x) if '(' not in x else x for x in mods ]
			if len(mods) != 0:
				return ' AND '.join(mods)

	if hasattr(builder, 'homology'):
		dct = { k.split('_mod_')[0]:v for k,v in builder.homology.org_cplx_homolog.items() if '_mod_' in k }
	else:
		dct = { v:k for k,v in builder.org.protein_mod[['Core_enzyme']].to_dict()['Core_enzyme'].items() }
	data['Cofactors in Modified Complex'] = data.apply(lambda x: cofactors(x, dct), axis = 1)

	# ME-model generics
	def generics_from_gene(x, dct):
		generics = []
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ]
		for key, lst in dct.items():
			for tag in [ x for y in tags for x in y ]:
				if tag + '-MONOMER' in lst:
					generics.append(key)
				if 'RNA_' + tag in lst:
					generics.append(key)

		if len(generics) != 0:
			return generics

	def generics_from_complex(x, dct):
		generics = []
		for key, lst in dct.items():
			if str(x['Complex ID']).split(':')[0] in lst:
				generics.append(key)
		if len(generics) != 0:
			return generics

	dct = { k.replace('generic_', ''):[ x.split('_mod_')[0] for x in v['enzymes'] ] for k,v in builder.org.generic_dict.items() }
	data['Generic Complex ID'] = data.apply(lambda x: generics_from_gene(x, dct), axis = 1)
	data['Generic Complex ID'].update(data.apply(lambda x: generics_from_complex(x, dct), axis = 1))
	data = data.explode('Generic Complex ID')

	# ME-model metacomplexes (e.g., ribosome)
	def rnapol(x, lst):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			if '{:s}-MONOMER'.format(tag) in lst:
				return 'RNAP'

	def ribosome(x, lst):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			if '{:s}-MONOMER'.format(tag) in lst or 'generic_{:s}'.format(tag) in lst:
				return 'ribosome:1'

	def degradosome(x, lst):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			if '{:s}-MONOMER'.format(tag) in lst:
				return 'RNA_degradosome:1'

	def excision(x, dct):
		subrxns = []
		for key, lst in dct.items():
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				if '{:s}-MONOMER'.format(tag) in lst or 'generic_{:s}'.format(tag) in lst:
					subrxns.append(key + ':1')
		if len(subrxns) != 0:
			return subrxns

	def transpaths(x, dct):
		pathways = []
		for key, lst in dct.items():
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				if '{:s}-MONOMER'.format(tag) in lst or 'generic_{:s}'.format(tag) in lst:
					pathways.append('translocation_pathway_' + key)
		if len(pathways) != 0:
			return pathways

	lst = builder.org.RNAP if isinstance(builder.org.RNAP, list) else [builder.org.RNAP]
	lst += list(builder.org.sigmas.index)
	data['MetaComplex ID'] = data.apply(lambda x: rnapol(x, lst), axis = 1)

	lst = [ x for y in [builder.org.ribosome_stoich[x]['stoich'].keys() for x in ['30_S_assembly', '50_S_assembly']] for x in y ]
	data['MetaComplex ID'].update(data.apply(lambda x: ribosome(x, lst), axis = 1))

	lst = [ x.split('_mod_')[0] for x in builder.org.rna_degradosome['rna_degradosome']['enzymes'] ]
	data['MetaComplex ID'].update(data.apply(lambda x: degradosome(x, lst), axis = 1))

	dct = { k:[ x.split('_mod_')[0] for x in v['enzymes'] ] for k,v in builder.org.excision_machinery.items() }
	data['MetaComplex ID'].update(data.apply(lambda x: excision(x, dct), axis = 1))

	dct = { k:list(v['enzymes'].keys()) for k,v in builder.org.translocation_pathways.items() if len(v['enzymes']) != 0 }
	data['MetaComplex ID'].update(data.apply(lambda x: transpaths(x, dct), axis = 1))
	data = data.explode('MetaComplex ID')

	# ME-model subreactions
	def ribosome_subrxns(x, dct):
		for key, lst in dct.items():
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				if '{:s}-MONOMER'.format(tag) in lst or 'generic_{:s}'.format(tag) in lst:
					return 'Ribosome_' + key

	def translation_subrxns(x, dct):
		for key, lst in dct.items():
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				if '{:s}-MONOMER'.format(tag) in lst or 'generic_{:s}'.format(tag) in lst:
					if key.endswith('InfA') or key.endswith('InfC'):
						return key
					elif key.endswith('InfB'):
						return 'Translation_initiation_gtp_factor_InfB'
					elif key == 'fmet_addition_at_START':
						return 'Translation_initiation_' + key
					elif key == 'FusA_mono_elongation':
						return 'Translation_elongation_FusA_mono'
					elif key == 'Tuf_gtp_regeneration':
						return 'Translation_elongation_' + key
					elif key in ['N_terminal_methionine_cleavage', 'DnaK_dependent_folding', 'GroEL_dependent_folding']:
						return 'Protein_processing_' + key
					elif key == 'PrfA_mono_mediated_termination':
						return 'Translation_termination_PrfA_mono_mediated'
					elif key == 'PrfB_mono_mediated_termination':
						return 'Translation_termination_PrfB_mono_mediated'
					elif key == 'generic_RF_mediated_termination':
						return 'Translation_termination_generic_RF_mediated'
					else:
						return 'Translation_termination_' + key

	def transcription_subrxns(x, dct):
		subrxns = []
		for key, lst in dct.items():
			tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
			tags = [ str(x).split(';') for x in tags ]
			for tag in [ x for y in tags for x in y ]:
				if '{:s}-MONOMER'.format(tag) in lst or 'generic_{:s}'.format(tag) in lst:
					subrxns.append(key)
		if len(subrxns) != 0:
			return subrxns

	dct = { k:[ x.split('_mod_')[0] for x in [v['enzyme']] ] for k,v in builder.org.ribosome_subreactions.items() }
	data['ME-model SubReaction'] = data.apply(lambda x: ribosome_subrxns(x, dct), axis = 1)

	dct = { k:[ x.split('_mod_')[0] for x in v['enzymes'] ] for k,v in builder.org.initiation_subreactions.items() }
	dct.update({ k:[ x.split('_mod_')[0] for x in v['enzymes'] ] for k,v in builder.org.elongation_subreactions.items() })
	dct.update({ k:[ x.split('_mod_')[0] for x in v['enzymes'] ] for k,v in builder.org.termination_subreactions.items() })
	data['ME-model SubReaction'].update(data.apply(lambda x: translation_subrxns(x, dct), axis = 1))

	dct = { k:[ x.split('_mod_')[0] for x in v['enzymes'] ] for k,v in builder.org.transcription_subreactions.items() }
	data['ME-model SubReaction'].update(data.apply(lambda x: transcription_subrxns(x, dct), axis = 1))
	data = data.explode('ME-model SubReaction')

	# Post-translation processing
	def processing(x, lst):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'] ]
		tags = [ str(x).split(';') for x in tags ]
		for tag in [ x for y in tags for x in y ]:
			if tag in lst:
				return 'TRUE'

	lst = builder.org.folding_dict['GroEL_dependent_folding']['enzymes']
	data['GroEL_dependent_folding'] = data.apply(lambda x: processing(x, lst), axis = 1)
	lst = builder.org.folding_dict['DnaK_dependent_folding']['enzymes']
	data['DnaK_dependent_folding'] = data.apply(lambda x: processing(x, lst), axis = 1)
	lst = builder.org.cleaved_methionine
	data['N_terminal_methionine_cleavage'] = data.apply(lambda x: processing(x, lst), axis = 1)

	# protein location
	def location(x, df):
		tags = [ x['Gene Locus ID'], x['Old Locus Tag'], x['BioCyc'], x['Generic Complex ID'] ]
		tags = [ str(x).split(';') for x in tags ] # we convert here None to 'None'
		tags = [ x for y in tags for x in y ]
		# extra modification if compared this function to others
		tags = [ '{:s}(\(\d*\)|\(\d*\:\d*\))'.format(x) for x in tags if x != 'None' ]

		#res = df[df['Protein'].str.match(x)][['Complex_compartment', 'Protein_compartment', 'translocase_pathway']].values
	    #query = [ '{:s}()'.format(y) for y in tags ] # sp old_locus_tag's, the query can be a string of tags separated by ';'
		res = df[df['Protein'].str.fullmatch('|'.join(tags))][['Complex_compartment', 'Protein_compartment', 'translocase_pathway']].values

		if len(res) != 0:
			return res[0]
		else:
			return [None, None, None]

	df = builder.org.protein_location.dropna(how = 'all', subset = ['Complex_compartment', 'Protein', 'Protein_compartment'])
	data['Complex Location'], data['Subunit Location'], data['Translocation Pathway'] = zip(*data.apply(lambda x: location(x, df), axis = 1))
	data = data.explode('Complex Location')

	# Filter in inferred enzyme-reaction associations
	cplxs = builder.org.enz_rxn_assoc_df.copy(deep = True)
	cplxs['Complexes'] = cplxs['Complexes'].apply(lambda gpr: [ x.split('_mod_')[0] for x in gpr.split(' OR ') ])
	cplxs = set([ '{:s}:\d+'.format(x) for x in cplxs.explode('Complexes')['Complexes'].to_list() ])
	# add complexes inferred from tRNA synthetases
	dct = builder.org.amino_acid_trna_synthetase
	cplxs.update([ '{:s}:\d+'.format(v.split('_mod_')[0]) for k,v in dct.items() ])

	# this dataframe contains only genes associated to a M-model reaction
	tmp1 = data.copy(deep = True).reset_index(drop = True)
	tmp1 = tmp1[tmp1['M-model Reaction ID'].notna() & tmp1['Complex ID'].notna()]
	tmp1 = tmp1[tmp1['Complex ID'].str.fullmatch('|'.join(cplxs))]

	# this dataframe contains genes associated to generics (correct association) and to reactions (incorrect association)
	tmp2 = data.copy(deep = True).reset_index(drop = True)
	tmp2 = tmp2[tmp2['M-model Reaction ID'].notna() & tmp2['Complex ID'].notna() & tmp2['Generic Complex ID'].notna()]
	tmp2 = tmp2[~tmp2['Complex ID'].str.fullmatch('|'.join(cplxs))]

	if not tmp2.empty:
		# TODO: Check stoichiometry of generics in complexes
		no_generics = tmp2.copy(deep = True)
		no_generics['Complex ID'] = no_generics.apply(lambda x: 'CPLX_{:s}-0:1({:s})'.format(x['M-model Reaction ID'], x['Generic Complex ID']), axis = 1)
		no_generics['Generic Complex ID'] = None
		no_generics.drop_duplicates(subset = ['Complex ID'])

		no_reactions = tmp2.copy(deep = True)
		no_reactions[['M-model Reaction ID', 'Reaction Name', 'Reversibility']] = None

		tmp2 = pandas.concat([no_reactions, no_generics], axis = 0)

	# this dataframe contains genes NOT associated to a M-model reaction
	tmp3 = data.copy(deep = True).reset_index(drop = True)
	tmp3 = tmp3[tmp3['M-model Reaction ID'].isna()]

	data = pandas.concat([tmp1, tmp2, tmp3], axis = 0)
	#data = data.drop_duplicates(inplace = False) # Is it correctly detecting duplicates when strings contain ()?

	def _save_to_excel(data, output):
		#if overwrite:
		try:
			pathlib.Path(output).unlink(missing_ok = True) # python>=3.8
		except:
			if pathlib.Path(output).exists():
				pathlib.Path(output).unlink() # python==3.7

		with open(output, 'wb') as outfile:
			writer = pandas.ExcelWriter(outfile, engine = 'xlsxwriter')
			data.to_excel(writer, index = False, freeze_panes = (1, 7))
			(max_row, max_col) = data.shape

			# Get the xlsxwriter workbook and worksheet objects.
			workbook  = writer.book
			worksheet = writer.sheets['Sheet1']

			# Set the autofilter.
			worksheet.autofilter(0, 0, max_row, max_col - 1)

			# Make the columns wider for clarity.
			worksheet.set_column_pixels(0,  max_col - 1, 96)

			# Close the Pandas Excel writer and output the Excel file.
			writer.close()
		return None

	try:
		# Save file as excel or tsv depending on the size
		if not output.endswith('.xlsx'):
			output = '.'.join(output.split('.')[:-1]) + '.xlsx'

		# GPRs can expand the model specifications beyond the max size of an excel file (1048576 rows)
		if data.shape[0] > 1048575: # one less to accommodate the header
			# Divide the DataFrame into pieces and save them
			rxn_ids = data.groupby(['M-model Reaction ID'], dropna = True).size() # Anything with 1 or more ocurrences

			gene_ids = data[data['M-model Reaction ID'].isin(rxn_ids[rxn_ids < 1000].index)]['Gene Locus ID']
			tmp = pandas.concat([ data[data['Gene Locus ID'].isin(gene_ids)], data[data['M-model Reaction ID'].isna()]], axis = 1)
			_save_to_excel(tmp, '.'.join(output.split('.')[:-1]) + '_{:02d}.xlsx'.format(0))

			if output.endswith('.xlsx'):
				for idx, gene_id in enumerate(rxn_ids[rxn_ids >= 1000].index):
					tmp = data[data['Gene Locus ID'].isin([gene_id])]
					_save_to_excel(tmp, '.'.join(output.split('.')[:-1]) + '_{:03d}.xlsx'.format(idx+1))

			# save the output as a tsv file
			if not output.endswith('.txt'):
				output = '.'.join(output.split('.')[:-2]) + '.txt'

			with open(output, 'w') as outfile:
				data.to_csv(outfile, index = False, sep = '\t')
		else:
			_save_to_excel(data, output)
	except:
		logging.warning('The builder.df_data was not saved to the \'{:s}\' file.'.format(output))

	return data

def correct_input(df):
	# correct Gene Locus ID to reflect if they are proteins or RNAs
	fn = lambda x: \
		'protein_{:s}'.format(x['Gene Locus ID']) if (x['Feature Type'] == 'CDS' and not x['Gene Locus ID'].startswith('protein_')) \
		else 'RNA_{:s}'.format(x['Gene Locus ID']) if (x['Feature Type'] != 'pseudo' and not x['Gene Locus ID'].startswith('RNA_')) \
		else 'pseudo_{:s}'.format(x['Gene Locus ID'])
	df['Gene Locus ID'] = df[['Gene Locus ID', 'Feature Type']].apply(fn, axis = 1)

	# correct generics to contain the `generic_` prefix
	fn = lambda x: 'generic_{:s}'.format(x) if isinstance(x, str) and not x.startswith('generic_') else x
	df['Generic Complex ID'] = df['Generic Complex ID'].apply(fn)
	return df

# 2nd step
def get_generics(df):
	#tmp = df[df['Generic Complex ID'].notna() & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['Generic Complex ID'].notna()]
	tmp = correct_input(tmp)

	fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
	tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)
	tmp = tmp.groupby(['Generic Complex ID']).agg({'Gene Locus ID': lambda x: sorted(set(x.tolist()))}).to_dict()['Gene Locus ID'].items()
	return tmp

def metacomplex_stoichiometry(df, key):
	#tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.startswith(key) & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.startswith(key)]
	if not tmp.empty:
		tmp = correct_input(tmp)

		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

		fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
			if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
		tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

		# collapse
		tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
		tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
		tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace

		tmp['stoich'] = tmp['MetaComplex ID'].str.split(':', expand = True).iloc[:, 1]
		return tmp
	else:
		return None

def ribosome_stoichiometry(df):
	tmp = metacomplex_stoichiometry(df, 'ribosome')
	ribosome_stoich = { 'ribosome' : { 'stoich' : { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich'])}}}
	ribosome_stoich['ribosome']['stoich']['gtp_c'] = 1
	return ribosome_stoich

def degradosome_stoichiometry(df):
	tmp = metacomplex_stoichiometry(df, 'RNA_degradosome')
	return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) } if tmp is not None else None

def dnapolymerase_stoichiometry(df):
	tmp = metacomplex_stoichiometry(df, 'DNAP')
	return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) }

def excision_machinery_stoichiometry(df, keys):
	#tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.match(keys) & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.match(keys)]
	if not tmp.empty:
		tmp = correct_input(tmp)

		fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
			if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
		tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

		# collapse
		tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
		tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
		tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace

		tmp['stoich'] = tmp['MetaComplex ID'].str.split(':', expand = True).iloc[:, 1]
		return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) }
	else:
		return None

def aa_synthetase_dict(df):
	#tmp = df[df['Definition'].str.contains('--tRNA ligase|-tRNA synthetase') & ~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['Definition'].str.contains('--tRNA ligase|-tRNA synthetase') & df['Definition'].notna()]
	tmp = correct_input(tmp)

	fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
	tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace

	tmp = tmp.drop_duplicates(['Gene Locus ID'])

	tmp = { k.split('--tRNA')[0].lower():v.split(':')[0] for k,v in zip(tmp['Definition'], tmp['Gene Locus ID']) }
	tmp = { k.replace('asparagine', 'asn'):v for k,v in tmp.items() }
	tmp = { k.replace('isoleucine', 'ile'):v for k,v in tmp.items() }
	tmp = { k.replace('tryptophan', 'trp'):v for k,v in tmp.items() }
	return { (k[:3]+'__L_c').replace('gly__L_c', 'gly_c'):v for k,v in tmp.items() }

def get_subreactions(df, key: str):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = correct_input(tmp)
	tmp = tmp[tmp['ME-model SubReaction'].notna()]

	tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}_cplx'.format(x))
	tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

	fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
	tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
	tmp = tmp.groupby(['ME-model SubReaction']).agg({'Gene Locus ID': lambda x: x.tolist()})
	return { k : {'enzymes' : list(set(v)), 'stoich': {}, 'element_contribution' : {}, 'keff' : []}
		for k,v in zip(tmp.index, tmp['Gene Locus ID']) if k.startswith(key) }

def get_df_rxns(df):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['M-model Reaction ID'].notna()]
	tmp = tmp[['M-model Reaction ID', 'Reaction Name', 'Reversibility']]
	tmp.columns = ['name', 'description', 'is_reversible']
	tmp = tmp.dropna().drop_duplicates('name', keep = 'first').set_index('name')

	return tmp

def get_df_cplxs(df, generics = False):
	tmp = df[df['Feature Type'].isin(['CDS', 'pseudo', 'ncRNA'])].fillna('')

	# get enzymatic complexes with generics subunits
	df_generics = df[df['Feature Type'].isin(['CDS'])].fillna('')
	df_generics = df_generics[df_generics['Complex ID'].str.contains('\(')]
	df_generics['Gene Locus ID'] = 'generic_' + df_generics['Complex ID'].str.extract('[A-Za-z0-9\_-]+:\d+\(([A-Za-z0-9\_-]+)\)', expand = True)
	df_generics['Complex ID'] = df_generics['Complex ID'].str.extract('([A-Za-z0-9\_-]+:\d+)\([A-Za-z0-9\_-]+\)', expand = True)

	tmp = pandas.concat([tmp, df_generics])
	tmp = tmp[~tmp['Complex ID'].str.contains('\(')]

	fn = lambda x: '{:s} complex'.format(x['Definition']) if x['Complex Name'] == '' else x['Complex Name']
	tmp['Complex Name'] = tmp[['Definition', 'Complex Name']].apply(fn, axis = 1)

	fn = lambda x: '{:s}({:s})'.format(x['Gene Locus ID'], x['Complex ID'].split(':')[1].split('(')[0]) if ':' in x['Complex ID'] else '{:s}(1)'.format(x['Gene Locus ID'])
	tmp['Gene Locus ID'] = tmp[['Gene Locus ID', 'Complex ID']].apply(fn, axis = 1)

	fn = lambda x: '{:s}-MONOMER'.format(x['Gene Locus ID'].split('(')[0]) if x['Complex ID'] == '' else x['Complex ID'].split('(')[0].split(':')[0]
	tmp['Complex ID'] = tmp[['Gene Locus ID', 'Complex ID']].apply(fn, axis = 1)

	fn = lambda x: x[0] if isinstance(x, list) else x.tolist()[0]
	tmp = tmp.groupby(['Complex ID']).agg({'Gene Locus ID': lambda x: ' AND '.join(sorted(set(x.tolist()))), 'Complex Name': fn}).reset_index()
	tmp = tmp[['Complex ID', 'Complex Name', 'Gene Locus ID']].set_index('Complex ID')

	return tmp

def get_df_ptms(df):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['Cofactors in Modified Complex'].notna()]

	if not tmp.empty:
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split('(')[0].split(':')[0] if isinstance(x, str) else x)

		fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')])
		tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

		tmp = tmp[['Modified Complex', 'Complex ID', 'Cofactors in Modified Complex']]
		tmp.columns = ['Modified_enzyme', 'Core_enzyme', 'Modifications']
		tmp = tmp.drop_duplicates(keep = 'first').set_index('Modified_enzyme')

		return tmp
	else:
		return pandas.DataFrame(columns = ['Modified Complex', 'Complex ID', 'Cofactors in Modified Complex'])

def get_df_enz2rxn(df, filter_in = set(), generics = False):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['M-model Reaction ID'].notna()]

	tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}-MONOMER'.format(x))
	tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

	fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace

	tmp = tmp.groupby(['M-model Reaction ID']).agg({'Gene Locus ID': lambda x: ' OR '.join(sorted(set(x.tolist())))})
	return tmp

def get_df_rna_enzs(df, filter_in = set(), generics = False):
	tmp = df[~df['Feature Type'].isin(['tRNA', 'rRNA'])]
	tmp = tmp[tmp['RNA mods/enzyme'].notna()]

	if not tmp.empty:
		tmp = correct_input(tmp)
		tmp['RNA mods/enzyme'] = tmp['RNA mods/enzyme'].str.split(',')
		tmp = tmp.explode('RNA mods/enzyme')

		tmp['modification'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[0])
		tmp['positions'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[1])

		tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}-MONOMER'.format(x))
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

		fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
			if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
		tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

		# collapse
		tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
		tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
		tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace

		tmp = tmp[['Gene Locus ID', 'modification', 'positions']]
		tmp.columns = ['enzymes', 'modification', 'positions']
		tmp = tmp.drop_duplicates(keep = 'first')

		return tmp
	else:
		return pandas.DataFrame(columns = ['enzymes', 'modification', 'positions'])

def get_df_rna_ptms(df, filter_in = set(), generics = False):
	tmp = df[df['Feature Type'].isin(['tRNA', 'rRNA'])]
	tmp = tmp[tmp['RNA mods/enzyme'].notna()]

	if not tmp.empty:
		tmp['RNA mods/enzyme'] = tmp['RNA mods/enzyme'].str.split(',')
		tmp = tmp.explode('RNA mods/enzyme')

		tmp['modification'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[0])
		tmp['positions'] = tmp['RNA mods/enzyme'].apply(lambda x: x.split('_at_')[1])

		# collapse
		tmp['Gene Locus ID'].update(tmp['Generic Complex ID']) # inplace

		tmp = tmp[['Gene Locus ID', 'modification', 'positions']]
		tmp.columns = ['bnum', 'modification', 'positions']
		tmp = tmp.drop_duplicates(keep = 'first')

		return tmp
	else:
		return pandas.DataFrame(columns = ['bnum', 'modification', 'positions'])

def get_df_protloc(df, filter_in = set(), generics = False):
	tmp = df.copy(deep = True)
	#tmp = tmp[~tmp['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['Complex Location'].notna()]

	if not tmp.empty:
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)
		tmp['Monomer ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}-MONOMER'.format(x))

		# collapse
		tmp['Monomer ID'].update(tmp['Complex ID']) # inplace

		tmp = tmp[['Monomer ID', 'Complex Location', 'Gene Locus ID', 'Subunit Location', 'Translocation Pathway']]
		tmp.columns = ['Complex', 'Complex_compartment', 'Protein', 'Protein_compartment', 'translocase_pathway']
		tmp = tmp.drop_duplicates(subset = ['Protein', 'Protein_compartment'], keep = 'first').set_index('Complex')

		return tmp
	else:
		return pandas.DataFrame(columns = ['Complex ID', 'Complex Location', 'Gene Locus ID', 'Subunit Location', 'Translocation Pathway'])

def get_df_transpaths(df, filter_in = set(), generics = False):
	tmp = df[~df['MetaComplex ID'].isnull()]
	tmp = tmp[tmp['MetaComplex ID'].str.match('translocation_pathway')]

	if not tmp.empty:
		tmp = correct_input(tmp)
		tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

		# collapse
		tmp['Complex ID'].update(tmp['Generic Complex ID'])
		tmp['Gene Locus ID'].update(tmp['Complex ID'])

		return tmp[['Gene Locus ID', 'MetaComplex ID']].groupby('MetaComplex ID').agg(lambda x: list(set(x.tolist())))
	else:
		return pandas.DataFrame(columns = ['Gene Locus ID'], index = ['MetaComplex ID'])

def get_df_input_from_excel(df, df_rxns):
	lst = [
		#'Complex Name',
		'Complex ID',
		'Cofactors in Modified Complex',
		'Generic Complex ID',
		'MetaComplex ID',
		'ME-model SubReaction',
		'M-model Reaction ID',
		'RNA mods/enzyme'
		]

	tmp1 = df[df['Feature Type'].str.match('CDS')].dropna(subset = lst, how = 'all', axis = 0)
	tmp2 = df[df['Feature Type'].isin(['rRNA', 'tRNA', 'ncRNA', 'tmRNA', 'misc_RNA', 'pseudo'])]
	df = pandas.concat([tmp1, tmp2], axis = 0).drop_duplicates()
	df = df.fillna({
		'Gene Locus ID' : '',
		'Reversibility' : False,
		#'Spontaneous' : False, # see reactions.txt input file
		'GroEL_dependent_folding' : False,
		'DnaK_dependent_folding' : False,
		'N_terminal_methionine_cleavage' : False
		})

	df_rxns = pandas.concat([get_df_rxns(df), df_rxns])
	df_cplxs = get_df_cplxs(df)
	df_ptms = get_df_ptms(df)
	df_enz2rxn = get_df_enz2rxn(df)
	df_protloc = get_df_protloc(df)
	df_transpaths = get_df_transpaths(df)

	# modifications of tRNAs and rRNAs
	df_rna_enzs = get_df_rna_enzs(df)
	df_rna_ptms = get_df_rna_ptms(df)
	cols = ['modification', 'positions']

	if not df_rna_ptms.empty or not df_rna_enzs.empty:
		df_rna_mods = pandas.merge(df_rna_ptms, df_rna_enzs, how = 'left', left_on = cols, right_on = cols).fillna('No_Machine')
	else:
		df_rna_mods = pandas.DataFrame(columns = ['bnum', 'modification', 'positions', 'enzymes'])

	return df, df_rxns, df_cplxs, df_ptms, df_enz2rxn, df_rna_mods, df_protloc, df_transpaths

def filter_out_cofactor(df, cofactors_to_skip):
	df = df.copy(deep = True)
	def filt(old_cofactors, *lst):
		if isinstance(old_cofactors, str):
			new_string = ' AND '.join([ cofactor for cofactor in old_cofactors.split(' AND ') if cofactor.split('(')[0] not in lst ])
			return new_string if new_string != '' else None
		else:
			return old_cofactors

	df['Cofactors in Modified Complex'] = df['Cofactors in Modified Complex'].apply(filt, args = (cofactors_to_skip))
	return df

def filter_in_cofactor(df, cofactors_to_keep):
	df = df.copy(deep = True)
	def filt(old_cofactors, *lst):
		if isinstance(old_cofactors, str):
			new_string = ' AND '.join([ cofactor for cofactor in old_cofactors.split(' AND ') if cofactor.split('(')[0] in lst ])
			return new_string if new_string != '' else None
		else:
			return old_cofactors

	df['Cofactors in Modified Complex'] = df['Cofactors in Modified Complex'].apply(filt, args = (cofactors_to_keep))
	return df
