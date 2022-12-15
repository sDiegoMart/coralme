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

def generate_organism_specific_matrix(genbank, m_model, output):
	contigs = []
	for contig in SeqIO.parse(genbank, 'genbank'):
		contigs.append(contig)

	# get all features
	lst = [ x for y in [ x.features for x in contigs ] for x in y ]
	lst = [ x for x in lst if x.type in [ 'CDS', 'ncRNA', 'tRNA', 'rRNA', 'tmRNA' ] ]

	# create a pandas DataFrame with organism-specific information to be completed with the builder data
	df = pandas.DataFrame(columns = [
		'Gene Locus ID',
		'Reference BBH',
		'Definition',
		'Feature Type',
		'Complex Name',
		'Complex ID',
		'Cofactors in Modified Complex',
		'Generic Complex ID',
		'MetaComplex ID',
		'ME-Model SubReaction',
		'M-Model Reaction ID',
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
		try:
			return [ x.id for x in m_model.genes.get_by_id(x).reactions ]
		except:
			return None

	def get_reaction_name(x):
		if x is not None:
			try:
				return m_model.reactions.get_by_id(x).name
			except:
				return None
		else:
			return None

	def get_reversibility(x):
		if x is not None:
			try:
				return 'TRUE' if m_model.reactions.get_by_id(x).reversibility else None
			except:
				return None
		else:
			return None

	def get_spontaneous(x):
		if x is not None:
			if 'spontaneous' in m_model.reactions.get_by_id(x).name:
				return 'TRUE'
			else:
				return None
		else:
			return None

	df['Gene Locus ID'] = [ x.qualifiers['locus_tag'][0] for x in lst ]
	df['Definition'] = [ x.qualifiers['product'][0] for x in lst ]
	df['Feature Type'] = [ x.type if x.qualifiers.get('pseudo') is None else 'pseudo' for x in lst ]
	df['M-Model Reaction ID'] = df['Gene Locus ID'].apply(get_reaction)
	df = df.explode('M-Model Reaction ID')
	df['Reaction Name'] = df['M-Model Reaction ID'].apply(get_reaction_name)
	df['Reversibility'] = df['M-Model Reaction ID'].apply(get_reversibility)
	# df['Spontaneous'] = df['M-Model Reaction ID'].apply(get_spontaneous)

	# df.set_index(['Gene Locus ID', 'Definition', 'Feature type'], inplace = True)
	df = df.sort_values(['M-Model Reaction ID', 'Gene Locus ID'])

	if pathlib.Path(output).is_file():
		pass
	else:
		df.to_excel(output, index = False)

	return df

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
	tmp = tmp[tmp['ME-Model SubReaction'].notna()]

	tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}_cplx'.format(x))
	tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

	fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	tmp['Modified Complex'].update(tmp['Generic Complex ID']) # inplace
	tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace
	tmp = tmp.groupby(['ME-Model SubReaction']).agg({'Gene Locus ID': lambda x: x.tolist()})
	return { k : {'enzymes' : list(set(v)), 'stoich': {}, 'element_contribution' : {}, 'keff' : []}
		for k,v in zip(tmp.index, tmp['Gene Locus ID']) if k.startswith(key) }

def get_df_rxns(df):
	tmp = df[~df['Feature Type'].isin(['pseudo'])]
	tmp = tmp[tmp['M-Model Reaction ID'].notna()]
	tmp = tmp[['M-Model Reaction ID', 'Reaction Name', 'Reversibility']]
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

	fn = lambda x: '{:s}_MONOMER'.format(x['Gene Locus ID'].split('(')[0]) if x['Complex ID'] == '' else x['Complex ID'].split('(')[0].split(':')[0]
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
	tmp = tmp[tmp['M-Model Reaction ID'].notna()]

	tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}_MONOMER'.format(x))
	tmp['Complex ID'] = tmp['Complex ID'].apply(lambda x: x.split(':')[0] if isinstance(x, str) else x)

	fn = lambda x: x['Complex ID'] + ''.join([ '_mod_{:s}'.format(x) for x in x['Cofactors in Modified Complex'].split(' AND ')]) \
		if isinstance(x['Cofactors in Modified Complex'], str) else numpy.nan
	tmp['Modified Complex'] = tmp[['Complex ID', 'Cofactors in Modified Complex']].apply(fn, axis = 1)

	# collapse
	tmp['Complex ID'].update(tmp['Modified Complex']) # inplace
	tmp['Gene Locus ID'].update(tmp['Complex ID']) # inplace

	tmp = tmp.groupby(['M-Model Reaction ID']).agg({'Gene Locus ID': lambda x: ' OR '.join(sorted(set(x.tolist())))})
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

		tmp['Gene Locus ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}_MONOMER'.format(x))
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
		tmp['Monomer ID'] = tmp['Gene Locus ID'].apply(lambda x: '{:s}_MONOMER'.format(x))

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
		'ME-Model SubReaction',
		'M-Model Reaction ID',
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
