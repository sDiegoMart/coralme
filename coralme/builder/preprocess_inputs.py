#!/usr/bin/python3
import numpy
import pathlib
import pandas
import warnings
from Bio import SeqIO

try:
	warnings.simplefilter(action = 'ignore', category = pandas.errors.SettingWithCopyWarning)
except:
	warnings.warn("This pandas version does not allow for correct warning handling. Pandas 1.5.1 is suggested.")

def generate_organism_specific_matrix(genbank, model):
	contigs = []
	for contig in SeqIO.parse(genbank, 'genbank'):
		contigs.append(contig)

	# get all features
	lst = [ x for y in [ x.features for x in contigs ] for x in y ]
	lst = [ x for x in lst if x.type in [ 'CDS', 'ncRNA', 'tRNA', 'rRNA', 'tmRNA' ] ]

	# create a pandas DataFrame with organism-specific information to be completed with the builder data
	df = pandas.DataFrame(columns = [
		'Gene Locus ID',
		'Gene Names',
		'Old Locus Tags',
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
		'Translocation Multiplier'
	])

	def get_reaction(x):
		if x is None:
			return None
		else:
			lst = []
			for gene in x.split(';'):
				if model.genes.has_id(gene):
					lst.append([ x.id for x in model.genes.get_by_id(gene).reactions ])

			# A (old)locustag/gene name can be associated to many reactions
			lst = [ x for y in lst for x in y ]
			if lst == []:
				return None
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

	df['Gene Locus ID'] = [ x.qualifiers.get('locus_tag', [None])[0] for x in lst ]
	df['Definition'] = [ x.qualifiers.get('product', [None])[0] for x in lst ]
	df['Feature Type'] = [ x.type if x.qualifiers.get('pseudo') is None else 'pseudo' for x in lst ]

	tmp = [ x.qualifiers.get('gene', None) for x in lst ]
	df['Gene Names'] = [ ';'.join(x) if x is not None else None for x in tmp ]

	tmp = [ x.qualifiers['old_locus_tag'] if x.qualifiers.get('old_locus_tag', None) is not None else None for x in lst ]
	df['Old Locus Tags'] = [ ';'.join(x) if x is not None else None for x in tmp ]

	df['M-model Reaction ID'] = df['Old Locus Tags'].apply(lambda x: get_reaction(x))
	df['M-model Reaction ID'].update(df['Gene Names'].apply(lambda x: get_reaction(x)))
	df['M-model Reaction ID'].update(df['Gene Locus ID'].apply(lambda x: get_reaction(x)))
	df = df.explode('M-model Reaction ID')

	df['Reaction Name'] = df['M-model Reaction ID'].apply(lambda x: get_reaction_name(x))
	df['Reversibility'] = df['M-model Reaction ID'].apply(lambda x: get_reversibility(x))

	# df.set_index(['Gene Locus ID', 'Definition', 'Feature type'], inplace = True)
	return df.sort_values(['M-model Reaction ID', 'Gene Locus ID'])

def complete_organism_specific_matrix(builder, data, model, output):
	def bbh(x, builder):
		if x is None:
			return None
		else:
			lst = []
			for gene in x.split(';'):
				lst.append(builder.homology.mutual_hits.get(gene, None))
			if set(lst) == set([None]):
				return None
			else:
				return [ x for x in lst if x is not None ][0]

	def monomers(x, builder):
		cplxID = builder.homology.org_cplx_homolog.get(x + '-MONOMER', None)
		if cplxID is None:
			return cplxID
		else:
			return cplxID + ':1'

	def complexes(x, df):
		cplxID = df[df['genes'].str.contains(x)]
		if len(cplxID) == 0:
			return None
		else:
			return ';'.join(cplxID.index.to_list())

	def cofactors(x, builder):
		dct = { k.split('_mod_')[0]:v for k,v in builder.homology.org_cplx_homolog.items() if '_mod_' in k }
		mods = dct.get(x + '-MONOMER', None)
		if mods is None:
			return mods
		else:
			return ' AND '.join([ x if '(' in x else x + '(1)' for x in mods.split('_mod_')[1:] ])

	def generics_from_gene(name, builder):
		lst = []
		for key, value in builder.org.generic_dict.items():
			values = [ x.split('_mod_')[0] for x in value['enzymes'] ]
			if name + '-MONOMER' in values:
				lst.append(key.replace('generic_', ''))
			if 'RNA_' + name in values:
				lst.append(key.replace('generic_', ''))

		if len(lst) >= 1:
			return lst

	def rnapol(x, builder):
		lst = [ x.replace('-MONOMER', '') for x in builder.org.RNAP ]
		if x['Gene Locus ID'] in lst:
			return 'RNAP'

	def ribosome(x, builder):
		lst = [ builder.org.ribosome_stoich[x]['stoich'].keys() for x in ['30_S_assembly', '50_S_assembly']]
		lst = [ x for y in lst for x in y ]

		if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
			return 'ribosome:1'

	def degradosome(x, builder):
		lst = builder.org.rna_degradosome['rna_degradosome']['enzymes']
		lst = [ x.split('_mod_')[0] for x in lst ]
		if x['Gene Locus ID'] + '-MONOMER' in lst:
			return 'RNA_degradosome:1'
		else:
			None

	def excision(x, builder):
		dct = builder.org.excision_machinery
		subrxns = []
		for key in dct.keys():
			lst = [ x.split('_mod_')[0] for x in dct[key]['enzymes'] ]
			if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
				subrxns.append(key + ':1')
		if len(subrxns) == 0:
			return None
		else:
			return subrxns

	def sigmas(x, builder):
		if isinstance(builder.org.RNAP, list):
			RNAP_components = builder.org.RNAP
		else:
			RNAP_components = [builder.org.RNAP]

		# combine
		lst = list(builder.org.sigmas.index) + RNAP_components
		if x['Gene Locus ID'] + '-MONOMER' in lst:
			return 'RNAP'

	def transpaths(x, builder):
		# simplified
		dct = { k:v['enzymes'] for k,v in builder.org.translocation_pathways.items() }
		pathways = []
		for key, value in dct.items():
			lst = value.keys()
			if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
				pathways.append('translocation_pathway_' + key)

		if len(pathways) != 0:
			return pathways
		else:
			return None

	def ribosome_subrxns(x, builder):
		for key, value in builder.org.ribosome_subreactions.items():
			lst = [value['enzyme']]
			lst = [ x.split('_mod_')[0] for x in lst ]
			if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
				return 'Ribosome_' + key

	def translation_subrxns(x, builder):
		for key, value in builder.org.initiation_subreactions.items():
			lst = value['enzymes']
			lst = [ x.split('_mod_')[0] for x in lst ]
			if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
				if 'InfA' in key or 'InfC' in key:
					return key
				elif key == 'Translation_gtp_initiation_factor_InfB':
					return 'Translation_initiation_gtp_factor_InfB'
				else:
					return 'Translation_initiation_' + key
		for key, value in builder.org.elongation_subreactions.items():
			lst = value['enzymes']
			lst = [ x.split('_mod_')[0] for x in lst ]
			if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
				if key == 'FusA_mono_elongation':
					return 'Translation_elongation_FusA_mono'
				else:
					return 'Translation_elongation_' + key
		for key, value in builder.org.termination_subreactions.items():
			lst = value['enzymes']
			lst = [ x.split('_mod_')[0] for x in lst ]
			if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
				if key in ['N_terminal_methionine_cleavage', 'DnaK_dependent_folding', 'GroEL_dependent_folding']:
					return 'Protein_processing_' + key
				elif key == 'PrfA_mono_mediated_termination':
					return 'Translation_termination_PrfA_mono_mediated'
				elif key == 'PrfB_mono_mediated_termination':
					return 'Translation_termination_PrfB_mono_mediated'
				elif key == 'generic_RF_mediated_termination':
					return 'Translation_termination_generic_RF_mediated'
				else:
					return 'Translation_termination_' + key

	def transcription_subrxns(x, builder):
		subrxns = []
		for key, value in builder.org.transcription_subreactions.items():
			lst = value['enzymes']
			lst = [ x.split('_mod_')[0] for x in lst ]
			if x['Gene Locus ID'] + '-MONOMER' in lst or 'generic_' + str(x['Generic Complex ID']) in lst:
				subrxns.append(key)
		if len(subrxns) == 0:
			return None
		else:
			return subrxns

	def groel(x, builder):
		if x in builder.org.folding_dict['GroEL_dependent_folding']['enzymes']:
			return 'TRUE'

	def dnak(x, builder):
		if x in builder.org.folding_dict['DnaK_dependent_folding']['enzymes']:
			return 'TRUE'

	def ntmeth(x, builder):
		if x in builder.org.cleaved_methionine:
			return 'TRUE'

	def location(x, builder):
		df = builder.org.protein_location.dropna(how = 'all', subset = ['Complex_compartment', 'Protein', 'Protein_compartment'])
		res = df[df['Protein'].str.match(x)][['Complex_compartment', 'Protein_compartment', 'translocase_pathway']].values
		if len(res) == 0:
			return (None, None, None)
		else:
			return res[0]

	def generics_from_complex(name, builder):
		if name is not None:
			name = name.split(':')[0]
			lst = []
			for key, value in builder.org.generic_dict.items():
				values = [ x.split('_mod_')[0] for x in value['enzymes'] ]

				if name in values:
					lst.append(key.replace('generic_', ''))

			if len(lst) >= 1:
				return lst

	def generics_in_complex(name, df):
		tmp = df[df['genes'].str.contains('generic')]
		if name is not None:
			lst = []
			if name in tmp.index:
				for generic in tmp['genes'][name]:
					lst.append('{:s}({:s})'.format(name, generic.replace('generic_', '').split('(')[0]))

			if len(lst) >= 1:
				return lst

	data['Reference BBH'] = data['Old Locus Tags'].apply(lambda x: bbh(x, builder))
	data['Reference BBH'].update(data['Gene Names'].apply(lambda x: bbh(x, builder)))
	data['Reference BBH'].update(data['Gene Locus ID'].apply(lambda x: bbh(x, builder)))

	dct = { v['Accession-1']:k for k,v in builder.org.gene_dictionary.iterrows() }

	def biocyc(x, builder):
		if x is None:
			return None
		else:
			lst = []
			for gene in x.split(';'):
				lst.append(dct.get(gene, None))
			if set(lst) == set([None]):
				return None
			else:
				return [ x for x in lst if x is not None ][0]

	data['BioCyc'] = data['Old Locus Tags'].apply(lambda x: biocyc(x, builder))
	data['BioCyc'].update(data['Gene Names'].apply(lambda x: biocyc(x, builder)))

	df = builder.org.complexes_df.copy(deep = True)
	#df = df[~df.index.str.contains('MONOMER')]
	df = df[df['genes'].str.contains('\(')]
	df['genes'] = df['genes'].str.split(' AND ')
	df = df.explode('genes')
	df['stoich'] = df['genes'].apply(lambda x: '1' if x.split('(')[1][:-1] == '' else str(x.split('(')[1][:-1]))
	df.index = df.index + ':' + df['stoich']

	# This overwrites the 'monomers' names with 'complexes' names
	#data['Complex ID'] = data['Gene Locus ID'].apply(lambda x: monomers(x, builder))
	#data['Complex ID'].update(data['Gene Locus ID'].apply(lambda x: complexes(x, df)))

	data['monomers'] = data['Gene Locus ID'].apply(lambda x: monomers(x, builder))
	data['complexes'] = data['Gene Locus ID'].apply(lambda x: complexes(x, df))

	def _combine(x):
		lst = [] if x['monomers'] is None else [x['monomers']]
		if x['complexes'] is None:
			pass
		else:
			for complex_name in x['complexes'].split(';'):
				lst.append(complex_name)
		return lst

	data['Complex ID'] = data[['monomers', 'complexes']].apply(lambda x: _combine(x), axis = 1)
	#data.drop('monomers', axis = 1, inplace = True)
	#data.drop('complexes', axis = 1, inplace = True)

	data['Cofactors in Modified Complex'] = data['Gene Locus ID'].apply(lambda x: cofactors(x, builder))
	data = data.explode('Complex ID')

	data['Generic Complex ID'] = data['Gene Locus ID'].apply(lambda x: generics_from_gene(x, builder))
	data = data.explode('Generic Complex ID')

	data['MetaComplex ID'] = data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: rnapol(x, builder), axis = 1)
	data['MetaComplex ID'].update(data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: ribosome(x, builder), axis = 1))
	data['MetaComplex ID'].update(data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: degradosome(x, builder), axis = 1))
	data['MetaComplex ID'].update(data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: excision(x, builder), axis = 1))
	data['MetaComplex ID'].update(data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: sigmas(x, builder), axis = 1))
	data['MetaComplex ID'].update(data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: transpaths(x, builder), axis = 1))
	data = data.explode('MetaComplex ID')

	data['ME-model SubReaction'] = data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: ribosome_subrxns(x, builder), axis = 1)
	data['ME-model SubReaction'].update(data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: translation_subrxns(x, builder), axis = 1))
	data['ME-model SubReaction'].update(data[['Gene Locus ID', 'Generic Complex ID']].apply(lambda x: transcription_subrxns(x, builder), axis = 1))
	data = data.explode('ME-model SubReaction')

	data['GroEL_dependent_folding'] = data['Gene Locus ID'].apply(lambda x: groel(x, builder))
	data['DnaK_dependent_folding'] = data['Gene Locus ID'].apply(lambda x: dnak(x, builder))
	data['N_terminal_methionine_cleavage'] = data['Gene Locus ID'].apply(lambda x: ntmeth(x, builder))

	data['Complex Location'], data['Subunit Location'], data['Translocation Pathway'] = zip(*data['Gene Locus ID'].apply(lambda x: location(x, builder)))
	data = data.explode('Complex Location')

	#data['Generic Complex ID'].update(data['Complex ID'].apply(lambda x: generics_from_complex(x, builder)))
	#data = data.explode('Generic Complex ID')

	#data['Complex ID'].update(data['Complex ID'].apply(lambda x: generics_in_complex(x, df)))
	#data = data.explode('Complex ID')

	if pathlib.Path(output).is_file():
		pass
	else:
		if output.endswith('xlsx'):
			pass
		else:
			output = '.'.join(output.split('.')[:-1]) + 'xlsx'

		with open(output, 'wb') as outfile:
			writer = pandas.ExcelWriter(outfile, engine = 'xlsxwriter')
			data.to_excel(writer, index = False, freeze_panes = (1, 2))
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
	tmp = df[df['Generic Complex ID'].notna() & ~df['Feature Type'].isin(['pseudo'])]
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
	tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.startswith(key) & ~df['Feature Type'].isin(['pseudo'])]
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

def ribosome_stoichiometry(df):
	tmp = metacomplex_stoichiometry(df, 'ribosome')
	ribosome_stoich = { 'ribosome' : { 'stoich' : { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich'])}}}
	ribosome_stoich['ribosome']['stoich']['gtp_c'] = 1
	return ribosome_stoich

def degradosome_stoichiometry(df):
	tmp = metacomplex_stoichiometry(df, 'RNA_degradosome')
	return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) }

def dnapolymerase_stoichiometry(df):
	tmp = metacomplex_stoichiometry(df, 'DNAP')
	return { k.split(':')[0]:int(v) for k,v in zip(tmp['Gene Locus ID'], tmp['stoich']) }

def excision_machinery_stoichiometry(df, keys):
	tmp = df[df['MetaComplex ID'].notna() & df['MetaComplex ID'].str.match(keys) & ~df['Feature Type'].isin(['pseudo'])]
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
	tmp = df[df['Definition'].str.contains('--tRNA ligase|-tRNA synthetase') & ~df['Feature Type'].isin(['pseudo'])]
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
	tmp = df[~df['Feature Type'].isin(['pseudo'])]
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
	tmp = df[~df['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['M-model Reaction ID'].notna()]
	tmp = tmp[['M-model Reaction ID', 'Reaction Name', 'Reversibility']]
	tmp.columns = ['name', 'description', 'is_reversible']
	tmp = tmp.dropna().drop_duplicates('name', keep = 'first').set_index('name')

	return tmp

def get_df_cplxs(df, generics = False):
	tmp = df[df['Feature Type'].isin(['CDS'])].fillna('')

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
	tmp = df[~df['Feature Type'].isin(['pseudo'])]
	tmp = df[df['Cofactors in Modified Complex'].notna()]

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
	tmp = df[~df['Feature Type'].isin(['pseudo'])]
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
	tmp = df[~df['Feature Type'].isin(['pseudo'])]
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
		'Complex Name',
		'Complex ID',
		'Cofactors in Modified Complex',
		'Generic Complex ID',
		'MetaComplex ID',
		'ME-model SubReaction',
		'M-model Reaction ID',
		'RNA mods/enzyme'
		]

	tmp1 = df[df['Feature Type'].str.match('CDS')].dropna(subset = lst, how = 'all', axis = 0)
	tmp2 = df[df['Feature Type'].isin(['rRNA', 'tRNA', 'ncRNA', 'tmRNA'])]
	df = pandas.concat([tmp1, tmp2], axis = 0).drop_duplicates()
	df = df.fillna({
		'Reversibility' : False,
		'Spontaneous' : False,
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
	if not df_rna_ptms.empty and not df_rna_enzs.empty:
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
