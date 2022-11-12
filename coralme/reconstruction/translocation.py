from collections import defaultdict
from coralme.core.processdata import PostTranslationData
from coralme.core.reaction import PostTranslationReaction

# Some proteins require different numbers of a complex in order to be
# translocated by a pathway
multipliers = {
	'YidC_MONOMER': {
		'b1855': 2., 'b3731': 2.
		},
	'TatE_MONOMER': {
		'b4072': 20.0, 'b1475': 5.166666666666667,
		'b1474': 5.166666666666667, 'b2817': 22.0,
		'b0973': 14.0, 'b0972': 14.0, 'b0997': 13.5,
		'b3893': 15.5, 'b0152': 21.0, 'b3894': 15.5,
		'b1589': 14.0, 'b1588': 15.0, 'b2938': 6.25,
		'b1587': 14.0, 'b1872': 26.0, 'b2997': 14.0,
		'b2996': 21.0, 'b2994': 14.0, 'b2435': 20.0,
		'b0894': 26.0, 'b0123': 23.0, 'b2206': 27.0,
		'b2205': 20.0
		},
	'TatA_MONOMER': {
		'b4072': 20.0, 'b1475': 5.166666666666667,
		'b1474': 5.166666666666667, 'b2817': 22.0,
		'b0973': 14.0, 'b0972': 14.0, 'b0997': 13.5,
		'b3893': 15.5, 'b0152': 21.0, 'b3894': 15.5,
		'b1589': 14.0, 'b1588': 15.0, 'b2938': 6.25,
		'b1587': 14.0, 'b1872': 26.0, 'b2997': 14.0,
		'b2996': 21.0, 'b2994': 14.0, 'b2435': 20.0,
		'b0894': 26.0, 'b0123': 23.0, 'b2206': 27.0,
		'b2205': 20.0
		},
	}

multipliers_protein_keys = defaultdict(dict)
for enzyme, value in multipliers.items():
	for bnum in value.keys():
		multipliers_protein_keys['protein_' + bnum][enzyme] = value[bnum]

mmol = 6.022e20  # number of molecules per mmol
nm2_per_m2 = 1e18  # used to convert nm^2 to m^2

# add_translocation_pathways helper function
def add_translocation_data_and_reaction(model, pathways, preprocessed_id, processed_id, compartment, peptide_data, membrane_constraints, alt = False):
	suffix = '_alt' if alt else ''

	data = PostTranslationData('translocation_' + preprocessed_id + '_' + compartment + suffix, model, processed_id, preprocessed_id)
	data.translocation = pathways
	data.translocation_multipliers = multipliers_protein_keys.get(preprocessed_id, {})

	# Add protein surface area constraint
	if membrane_constraints and compartment != 'Periplasm':
		protein = peptide_data.protein
		protein_met = model.metabolites.get_by_id('protein_' + protein)
		mass = protein_met.formula_weight / 1000.  # in kDa
		membrane_thickness = model.global_info['membrane_thickness']
		thickness = membrane_thickness[compartment]
		# Relationship uses protein molecular in kDa
		# Adds surface area constraint in units of m^2/mmol
		data.surface_area['SA_protein_' + compartment] = (1.21 / thickness * 2.) * mass * mmol / nm2_per_m2

	rxn = PostTranslationReaction('translocation_' + peptide_data.id + '_' + compartment + suffix)
	rxn.posttranslation_data = data
	model.add_reaction(rxn)
	rxn.update()

def add_translocation_pathways(model, pathways_df, abbreviation_to_pathway, membrane_constraints = False):
	# loop through all translation data and add translocation rxns/surface area constraints if they are membrane proteins
	# We can save time here if filtering what we need to process, not iterate over all TranslationData
	#for peptide_data in model.translation_data:
	for idx, translocation_info in pathways_df.iterrows():
		peptide_data = model.process_data.get_by_id(translocation_info['Protein'])
		# extract translocation info if peptide contained in complex stoichiometry
		#translocation_info = pathways_df[pathways_df['Protein'].str.match(peptide_data.id)]
		# iterate if protein is not in a membrane complex
		#if len(translocation_info) == 0:
			#continue

		# Assign preprocessed and processed (translocated) peptide ids
		compartment = translocation_info['Protein_compartment'] #.values[0]
		processed_id = 'protein_' + peptide_data.id + '_' + compartment
		preprocessed_id = 'protein_' + peptide_data.id

		# compile translocation pathways for each membrane protein
		pathways = set()
		pathways_alt = set()
		for abbrev in translocation_info['translocase_pathway']: #.values[0]:
			try:
				pathway_name = abbreviation_to_pathway[abbrev]
			except:
				continue

			# The tat translocation pathway can use an alternate enzyme
			if isinstance(pathway_name, list):
				pathways.add(pathway_name[0])
				pathways_alt.add(pathway_name[1])
			else:
				pathways.add(pathway_name)
				pathways_alt.add(pathway_name)

		# Call the helper function
		add_translocation_data_and_reaction(
			model, pathways, preprocessed_id, processed_id, compartment, peptide_data, membrane_constraints, alt = False)
		# if there's an alternative pathway (tat) add this reaction as well
		if pathways != pathways_alt:
			add_translocation_data_and_reaction(
				model, pathways_alt, preprocessed_id, processed_id, compartment, peptide_data, membrane_constraints, alt = True)

# add_lipoprotein_formation helper function
def add_lipoprotein_data_and_reaction(first_lipid, second_lipid, update):

	# Add PostTranslation Data, modifications and surface area
	data = PostTranslationData(reaction_prefix + '_' + second_lipid, model, processed_id, preprocessed_id)
	data.subreactions['mod_' + first_lipid] = 1
	data.subreactions['mod2_' + second_lipid + '_p'] = 1
	data.biomass_type = 'lipid_biomass'

	if membrane_constraints:
		thickness_dict = model.global_info['membrane_thickness']
		thickness = thickness_dict['Outer_Membrane']

		# From Liu et al. x2 for each to account for each leaflet
		protein_SA = 1.21 / thickness * 2 * mass * mmol / nm2_per_m2
		data.surface_area = {
			'SA_protein_' + compartment: -protein_SA,
			'SA_lipoprotein': 1. * mmol / nm2_per_m2
			}

	# Add Reaction to model and associated it with its data
	rxn = PostTranslationReaction(reaction_prefix + '_' + second_lipid)
	model.add_reaction(rxn)
	rxn.posttranslation_data = data

	if update:
		rxn.update()

def add_lipoprotein_formation(model, compartment_dict, lipoprotein_precursors, lipid_modifications, membrane_constraints = False, update = True):
	# loop through all proteins which need lipid modifications (lipoproteins)
	for protein in lipoprotein_precursors.values():
		compartment = compartment_dict.get(protein)
		protein_met = model.metabolites.get_by_id('protein_' + protein)
		mass = protein_met.formula_weight / 1000.  # in kDa

		processed_id = 'protein_' + protein + '_lipoprotein_' + compartment
		preprocessed_id = 'protein_' + protein + '_' + compartment

		for mod in lipid_modifications:
			reaction_prefix = protein + '_lipid_modification_' + mod
			add_lipoprotein_data_and_reaction(mod, 'pg160', update = update)
			add_lipoprotein_data_and_reaction(mod, 'pe160', update = update)
