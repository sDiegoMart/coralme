import re
import typing
import logging

# install by the user
import tqdm
bar_format = '{desc:<75}: {percentage:.1f}%|{bar}| {n_fmt:>5}/{total_fmt:>5} [{elapsed}<{remaining}]'
import numpy
import scipy
import sympy
import cobra
from cobra.medium import find_boundary_types, find_external_compartment, sbo_terms

import coralme

def _update(MEReaction):
	"""updates all component reactions"""
	MEReaction.update()
	return None

class MEModel(cobra.core.model.Model):
	def __init__(self, name, mu = 'mu'):
		cobra.Model.__init__(self, name)

		self.global_info = {
			'kt' : 4.5,
			'r0' : 0.087,
			'k_deg' : 12,
			'm_rr' : 1453.0,
			'm_aa' : 0.109,
			'm_nt' : 0.324,
			'f_rRNA' : 0.86,
			'f_mRNA' : 0.02,
			'f_tRNA' : 0.12,
			'm_tRNA' : 25.000,
			'temperature' : 37,
			'propensity_scaling' : 0.45,

			'dnapol_id' : 'DNAP',
			'ribosome_id' : 'ribosome',
			'dummy_rxn_id' : 'dummy_reaction',
			'degradosome_id' : 'RNA_degradosome',
			'mg2_per_ribosome' : 171,
			'amino_acid_loader' : 'generic_Tuf',
			'feature_types' : [ 'CDS', 'rRNA', 'tRNA', 'ncRNA', 'tmRNA' ],

			'translation_stop_dict' : {
				'UAG': 'PrfA_mono',
				'UGA': 'PrfB_mono',
				'UAA': 'generic_RF',
				},

			'transcription_subreactions' : {
				'Transcription_normal_rho_independent' : '',
				'Transcription_normal_rho_dependent' : 'atp_hydrolysis_rho',
				'Transcription_stable_rho_independent' : '',
				'Transcription_stable_rho_dependent' : 'atp_hydrolysis_rho',
				},

			'translation_subreactions' : {
				'Translation_initiation_factor_InfA' : '',
				'Translation_initiation_factor_InfC' : '',
				'Translation_initiation_fmet_addition_at_START' : 'FMETTRS',
				'Translation_initiation_gtp_factor_InfB' : 'atp_hydrolysis',
				'Translation_elongation_FusA_mono' : 'atp_hydrolysis',
				'Translation_elongation_Tuf_gtp_regeneration' : '',
				'Translation_termination_PrfA_mono_mediated' : '',
				'Translation_termination_PrfB_mono_mediated' : '',
				'Translation_termination_generic_RF_mediated' : '',
				'Translation_termination_peptide_deformylase_processing' : 'DEF',
				'Translation_termination_peptide_chain_release' : 'gtp_hydrolysis',
				'Translation_termination_ribosome_recycler' : '',
				'Protein_processing_GroEL_dependent_folding' : 'atp_hydrolysis_groel',
				'Protein_processing_DnaK_dependent_folding' : 'atp_hydrolysis',
				'Protein_processing_N_terminal_methionine_cleavage' : 'MAP',
				'Ribosome_RbfA_mono_assembly_factor_phase1' : '',
				'Ribosome_RimM_mono_assembly_factor_phase1' : '',
				'Ribosome_gtp_bound_30S_assembly_factor_phase1' : 'gtp_hydrolysis_era'
				},

			'peptide_processing_subreactions' : [
				'Translation_termination_peptide_chain_release',
				'Translation_termination_peptide_deformylase_processing',
				'Translation_termination_ribosome_recycler'
				],

			'translocation_pathway' : {
				'sec' : {
					'abbrev' : 's',
					'keff' : 4.0000,
					'length_dependent_energy' : 'True',
					'stoichiometry' : 'atp_hydrolysis_sec_pathway'
					},
				'secA' : {
					'abbrev' : 'a',
					'keff' : 4.0000,
					'length_dependent_energy' : 'True',
					'stoichiometry' : 'atp_hydrolysis_secA'
					},
				'tat' : {
					'abbrev' : 't',
					'keff' : 0.0125,
					'length_dependent_energy' :
					'False', 'stoichiometry' : ''
					},
				'tat_alt' : {
					'abbrev' : 't',
					'keff' : 0.0125,
					'length_dependent_energy' :
					'False', 'stoichiometry' : ''
					},
				'yidC' : {
					'abbrev' : 'y',
					'keff' : 20.000,
					'length_dependent_energy' : 'False',
					'stoichiometry' : 'gtp_hydrolysis'
					},
				'srp' : {
					'abbrev' : 'r',
					'keff' : 20.000,
					'length_dependent_energy' : 'False',
					'stoichiometry' : 'gtp_hydrolysis_srp_pathway',
					'FtsY' : 'FtsY_MONOMER'
					},
				'srp_yidC' : {
					'abbrev' : 'p',
					'keff' : 20.000,
					'length_dependent_energy' : 'False',
					'stoichiometry' : 'gtp_hydrolysis'
					},
				'lol' : {
					'abbrev' : 'l',
					'keff' : 0.9000,
					'length_dependent_energy' : 'False',
					'stoichiometry' : 'atp_hydrolysis'
					},
				'bam' : {
					'abbrev' : 'b',
					'keff' : 0.0270,
					'length_dependent_energy' : 'False',
					'stoichiometry' : ''
					}
				},

			'excision_machinery' : [
				'rRNA_containing',
				'monocistronic',
				'polycistronic_wout_rRNA'
				],

			'biomass_constraints' : [
				'protein_biomass',
				'mRNA_biomass',
				'tRNA_biomass',
				'rRNA_biomass',
				'ncRNA_biomass',
				'tmRNA_biomass',
				'DNA_biomass',
				'lipid_biomass',
				'constituent_biomass',
				'prosthetic_group_biomass',
				'peptidoglycan_biomass'
				],

			'compartments' : {
				'c'    : 'Cytoplasm',
				'e'    : 'Extracellular',
				'p'    : 'Periplasm',
				'mc'   : 'ME-Model Constraint'
				},

			'START_tRNA' : [],
			'rna_components' : [],
			'knockouts' : [],
			'genome_mods' : {},
			'trna_misacylation' : {},

			'me.gam' : 45.,
			'me.ngam' : 1.,
			'me.unmodeled_protein_fraction' : 0.36
			}

		self.process_data = cobra.core.dictlist.DictList()
		self.metabolites = cobra.core.dictlist.DictList()

		# set growth rate symbolic variable
		self.mu = sympy.Symbol(mu, positive = True)

		# Create the biomass dilution constraint
		self._biomass = coralme.core.component.Constraint('biomass')
		self._biomass_dilution = coralme.core.reaction.SummaryVariable('biomass_dilution')
		self._biomass_dilution.add_metabolites({self._biomass: -1})
		self.add_reactions([self._biomass_dilution])

		# cobra/core/reaction.py:328 Cannot convert expression to float
		# Override method: solved to check if variable type is sympy.core.symbol.Symbol or float
		self._biomass_dilution.upper_bound = self.mu
		self._biomass_dilution.lower_bound = self.mu

		# Maintenance energy
		self._gam = self.global_info['me.gam'] # default value
		self._ngam = self.global_info['me.ngam'] # default value

		"""
		Unmodeled protein is handled by converting protein_biomass to
		biomass, and requiring production of the appropriate amount of dummy
		protein
		"""
		self._unmodeled_protein_fraction = self.global_info['me.unmodeled_protein_fraction'] # default value

	def add_boundary(
		self,
		metabolite,
		type="exchange",
		reaction_id=None,
		lb=None,
		ub=None,
		sbo_term=None,
	):
		"""
		Add a boundary reaction for a given metabolite.

		There are three different types of pre-defined boundary reactions:
		exchange, demand, and sink reactions.
		An exchange reaction is a reversible, unbalanced reaction that adds
		to or removes an extracellular metabolite from the extracellular
		compartment.
		A demand reaction is an irreversible reaction that consumes an
		intracellular metabolite.
		A sink is similar to an exchange but specifically for intracellular
		metabolites, i.e., a reversible reaction that adds or removes an
		intracellular metabolite.

		If you set the reaction `type` to something else, you must specify the
		desired identifier of the created reaction along with its upper and
		lower bound. The name will be given by the metabolite name and the
		given `type`.

		Parameters
		----------
		metabolite : cobra.Metabolite
			Any given metabolite. The compartment is not checked but you are
			encouraged to stick to the definition of exchanges and sinks.
		type : str, {"exchange", "demand", "sink"}
			Using one of the pre-defined reaction types is easiest. If you
			want to create your own kind of boundary reaction choose
			any other string, e.g., 'my-boundary'.
		reaction_id : str, optional
			The ID of the resulting reaction. This takes precedence over the
			auto-generated identifiers but beware that it might make boundary
			reactions harder to identify afterwards when using `model.boundary`
			or specifically `model.exchanges` etc.
		lb : float, optional
			The lower bound of the resulting reaction.
		ub : float, optional
			The upper bound of the resulting reaction.
		sbo_term : str, optional
			A correct SBO term is set for the available types. If a custom
			type is chosen, a suitable SBO term should also be set.

		Returns
		-------
		cobra.Reaction
			The created boundary reaction.

		Examples
		--------
		>>> from cobra.io load_model
		>>> model = load_model("textbook")
		>>> demand = model.add_boundary(model.metabolites.atp_c, type="demand")
		>>> demand.id
		'DM_atp_c'
		>>> demand.name
		'ATP demand'
		>>> demand.bounds
		(0, 1000.0)
		>>> demand.build_reaction_string()
		'atp_c --> '

		"""
		ub = +1000 if ub is None else ub
		lb = -1000 if lb is None else lb
		types = {
			"exchange": ("EX", lb, ub, sbo_terms["exchange"]),
			"demand": ("DM", 0, ub, sbo_terms["demand"]),
			"sink": ("SK", lb, ub, sbo_terms["sink"]),
		}
		if type == "exchange":
			external = find_external_compartment(self)
			if metabolite.compartment != external:
				raise ValueError(
					"The metabolite is not an external metabolite"
					" (compartment is `%s` but should be `%s`). "
					"Did you mean to add a demand or sink? "
					"If not, either change its compartment or "
					"rename the model compartments to fix this."
					% (metabolite.compartment, external)
				)
		if type in types:
			prefix, lb, ub, default_term = types[type]
			if reaction_id is None:
				reaction_id = "{}_{}".format(prefix, metabolite.id)
			if sbo_term is None:
				sbo_term = default_term
		if reaction_id is None:
			raise ValueError(
				"Custom types of boundary reactions require a custom "
				"identifier. Please set the `reaction_id`."
			)
		if reaction_id in self.reactions:
			raise ValueError(
				"Boundary reaction '{}' already exists.".format(reaction_id)
			)
		name = "{} {}".format(metabolite.name, type)
		rxn = cobra.core.reaction.Reaction(id=reaction_id, name=name, lower_bound=lb, upper_bound=ub)
		rxn.add_metabolites({metabolite: -1})
		if sbo_term:
			rxn.annotation["sbo"] = sbo_term
		cobra.core.model.Model.add_reactions(self, [rxn])
		return rxn

	def add_reactions(self, reaction_list):
		"""Add reactions to the model.

		Reactions with identifiers identical to a reaction already in the
		model are ignored.

		This method was modified from the original cobrapy to not populate
		the solver interface, effectively speeding up the reconstruction.

		Parameters
		----------
		reaction_list : list
			A list of `cobra.Reaction` objects
		"""

		def existing_filter(rxn):
			if rxn.id in self.reactions:
				return False
			return True

		# First check whether the reactions exist in the model.
		pruned = cobra.core.dictlist.DictList(filter(existing_filter, reaction_list))

		# Add reactions. Also take care of genes and metabolites in the loop.
		for reaction in pruned:
			reaction._model = self

			# Build a `list()` because the dict will be modified in the loop.
			for metabolite in list(reaction.metabolites):
				# TODO: Should we add a copy of the metabolite instead?
				if metabolite not in self.metabolites:
					self.add_metabolites(metabolite)
				# A copy of the metabolite exists in the model, the reaction
				# needs to point to the metabolite in the model.
				else:
					# FIXME: Modifying 'private' attributes is horrible.
					stoichiometry = reaction._metabolites.pop(metabolite)
					model_metabolite = self.metabolites.get_by_id(metabolite.id)
					reaction._metabolites[model_metabolite] = stoichiometry
					model_metabolite._reaction.add(reaction)

			for gene in list(reaction._genes):
				# If the gene is not in the model, add it
				if not self.genes.has_id(gene.id):
					self.genes += [gene]
					gene._model = self
				# Otherwise, make the gene point to the one in the model
				else:
					model_gene = self.genes.get_by_id(gene.id)
					if model_gene is not gene:
						reaction._dissociate_gene(gene)
						reaction._associate_gene(model_gene)

		self.reactions += pruned

		# from cameo ...
		#self._populate_solver(pruned)

	def add_biomass_constraints_to_model(self, biomass_types):
		for biomass_type in tqdm.tqdm(biomass_types, 'Adding biomass constraint(s) into the ME-Model...', bar_format = bar_format):
			if '_biomass' not in biomass_type:
				raise ValueError('Biomass types should be suffixed with \'_biomass\'.')
			constraint_obj = coralme.core.component.Constraint(biomass_type)
			summary_variable_obj = coralme.core.reaction.SummaryVariable('{:s}_to_biomass'.format(biomass_type))
			summary_variable_obj.add_metabolites({constraint_obj: -1, self._biomass: 1})
			self.add_reactions([summary_variable_obj])

	@property
	def unmodeled_protein(self):
		return self.metabolites.get_by_id('protein_dummy')

	@property
	def unmodeled_protein_biomass(self):
		return self.metabolites.get_by_id('unmodeled_protein_biomass')

	@property
	def unmodeled_protein_fraction(self):
		return self._unmodeled_protein_fraction

	@unmodeled_protein_fraction.setter
	def unmodeled_protein_fraction(self, value):
		if 'protein_biomass_to_biomass' not in self.reactions:
			raise UserWarning(
				'Must add SummaryVariable handling the protein '
				'biomass constraint (via :meth:`add_biomass_constraints_to_model`) '
				'before defining the unmodeled protein fraction'
				)

		# See the Biomass_formulations for an explanation
		amount = value / (1 - value)
		self.reactions.protein_biomass_to_biomass.add_metabolites({self.unmodeled_protein_biomass: -amount}, combine = False)
		self.reactions.protein_biomass_to_biomass.add_metabolites({self._biomass: 1 + amount}, combine = False)
		self._unmodeled_protein_fraction = value

	@property
	def gam(self):
		return self._gam

	@gam.setter
	def gam(self, value):
		if 'GAM' not in self.reactions:
			logging.warning('Adding GAM (ATP requirement for growth) reaction into the ME-Model.')
			self.add_reactions([coralme.core.reaction.SummaryVariable('GAM')])
			self.reactions.GAM.lower_bound = self.mu
		#atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1, 'pi_c': 1} # charges: -4, 0 => -3, +1, -2
		atp_hydrolysis = self.process_data.get_by_id('atp_hydrolysis').stoichiometry
		for met, coeff in atp_hydrolysis.items():
			self.reactions.GAM.add_metabolites({met: value * coeff}, combine = False)
		self._gam = value

	@property
	def ngam(self):
		return self._ngam

	@ngam.setter
	def ngam(self, value):
		if 'ATPM' not in self.reactions:
			logging.warning('Adding ATPM (ATP requirement for maintenance) reaction into the ME-Model.')
			#atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1, 'pi_c': 1} # charges: -4, 0 => -3, +1, -2
			atp_hydrolysis = self.process_data.get_by_id('atp_hydrolysis').stoichiometry
			self.add_reactions([coralme.core.reaction.SummaryVariable('ATPM')])
			self.reactions.ATPM.add_metabolites(atp_hydrolysis)
		self.reactions.ATPM.lower_bound = value
		self._ngam = value

	# data types generators:
	# StoichiometricData, ComplexData, TranslationData, TranscriptionData,
	# GenericData, tRNAData, TranslocationData, PostTranslationData, SubreactionData
	@property
	def stoichiometric_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.StoichiometricData):
				yield data

	@property
	def complex_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.ComplexData):
				yield data

	@property
	def translation_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.TranslationData):
				yield data

	@property
	def transcription_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.TranscriptionData):
				yield data

	@property
	def generic_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.GenericData):
				yield data

	@property
	def tRNA_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.tRNAData):
				yield data

	@property
	def translocation_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.TranslocationData):
				yield data

	@property
	def posttranslation_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.PostTranslationData):
				yield data

	@property
	def subreaction_data(self):
		for data in self.process_data:
			if isinstance(data, coralme.core.processdata.SubreactionData):
				yield data

	# MEModel methods
	def get_metabolic_flux(self, solution = None):
		"""Extract the flux state for Metabolic reactions."""
		if solution is None:
			solution = self.solution
		if solution.status != 'optimal':
			raise ValueError('Solution status \'{:s}\' is not \'optimal\'.'.format(solution.status))
		flux_dict = {r.id: 0 for r in tqdm.tqdm(list(self.stoichiometric_data), 'Building reaction dictionary...', bar_format = bar_format)}
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-Model Reactions...', bar_format = bar_format):
			if isinstance(reaction, coralme.core.reaction.MetabolicReaction):
				m_reaction_id = reaction.stoichiometric_data.id
				if reaction.reverse:
					flux_dict[m_reaction_id] -= solution.fluxes[reaction.id]
				else:
					flux_dict[m_reaction_id] += solution.fluxes[reaction.id]
			# SummaryVariable in M-Model
			elif reaction.id == 'ATPM':
				flux_dict[reaction.id] = solution.fluxes[reaction.id]
			# Exchange, Demand, and Sink reactions
			elif reaction.id.startswith('EX_') or reaction.id.startswith('DM_') or reaction.id.startswith('SK_'):
				flux_dict[reaction.id] = solution.fluxes[reaction.id]
		return flux_dict

	def get_transcription_flux(self, solution = None):
		"""Extract the flux state of Transcription reactions."""
		if solution is None:
			solution = self.solution
		if solution.status != 'optimal':
			raise ValueError('Solution status \'{:s}\' is not \'optimal\'.'.format(solution.status))
		flux_dict = {}
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-Model Reactions...', bar_format = bar_format):
			if isinstance(reaction, coralme.core.reaction.TranscriptionReaction):
				for rna_id in reaction.transcription_data.RNA_products:
					locus_id = rna_id.replace('RNA_', '', 1)
					if locus_id not in flux_dict:
						flux_dict[locus_id] = 0
					flux_dict[locus_id] += solution.fluxes[reaction.id]
		return flux_dict

	def get_translation_flux(self, solution = None):
		"""Extract the flux state of Translation reactions."""
		if solution is None:
			solution = self.solution
		if solution.status != 'optimal':
			raise ValueError('Solution status \'{:s}\' is not \'optimal\'.'.format(solution.status))
		flux_dict = {r.id: 0 for r in tqdm.tqdm(list(self.translation_data), 'Building reaction dictionary...', bar_format = bar_format)}
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-Model Reactions...', bar_format = bar_format):
			if isinstance(reaction, coralme.core.reaction.TranslationReaction):
				protein_id = reaction.translation_data.id
				flux_dict[protein_id] += solution.fluxes[reaction.id]
		return flux_dict

	def construct_s_matrix(self, growth_rate) -> scipy.sparse.dok_matrix:
		"""Build the stoichiometric matrix at a specific growth rate."""
		# initialize to 0
		s_matrix = scipy.sparse.dok_matrix((len(self.metabolites), len(self.reactions)))
		# populate with stoichiometry
		for idx, rxn in tqdm.tqdm(enumerate(self.reactions), bar_format = bar_format):
			for met, value in rxn._metabolites.items():
				met_index = self.metabolites.index(met)
				if hasattr(value, 'subs'):
					s_matrix[met_index, idx] = float(value.subs(self.mu, growth_rate))
				else:
					s_matrix[met_index, idx] = float(value)
		return s_matrix

	def _construct_attribute_vector(self, attr_name, growth_rate):
		"""
		Build a vector of a reaction attribute at a specific growth rate.
		Mainly used for upper and lower bounds.
		"""
		return numpy.array([
			float(value.subs(self.mu, growth_rate))
			if hasattr(value, 'subs') else float(value)
			for value in tqdm.tqdm(self.reactions.list_attr(attr_name), bar_format = bar_format)
			])

	def compute_solution_error(self, solution = None):
		errors = {}
		if solution is None:
			solution = self.solution
		s_matrix = self.construct_s_matrix(solution.f)
		lb = self._construct_attribute_vector('lower_bound', solution.f)
		ub = self._construct_attribute_vector('upper_bound', solution.f)
		# old code
		#x = numpy.array(solution.x)
		x = numpy.array(list(solution.fluxes.values()))
		err = abs(s_matrix * x)
		errors['max_error'] = err.max()
		errors['sum_error'] = err.sum()
		ub_err = min(ub - x)
		errors['upper_bound_error'] = abs(ub_err) if ub_err < 0 else 0
		lb_err = min(x - lb)
		errors['lower_bound_error'] = abs(lb_err) if lb_err < 0 else 0
		return errors

	def prune(self, skip = None):
		"""
		Remove all unused metabolites and reactions
		This should be run after the model is fully built. It will be
		difficult to add new content to the model once this has been run.
		skip: list
			List of complexes/proteins/mRNAs/TUs to remain unpruned from model.
		"""
		if not skip:
			skip = []

		complex_data_list = [ i.id for i in self.complex_data if i.id not in skip ]

		for c_d in tqdm.tqdm(complex_data_list, 'Pruning unnecessary ComplexData reactions...', bar_format = bar_format):
			c = self.process_data.get_by_id(c_d)
			cplx = c.complex
			if len(cplx.reactions) == 1:
				list(cplx.reactions)[0].delete(remove_orphans = True)
				self.process_data.remove(self.process_data.get_by_id(c_d))

		for p in tqdm.tqdm(list(self.metabolites.query('_folded')), 'Pruning unnecessary FoldedProtein reactions...', bar_format = bar_format):
			if 'partially' not in p.id and p.id not in skip:
				delete = True
				for rxn in p.reactions:
					if rxn.metabolites[p] < 0:
						delete = False
						break

				if delete:
					while len(p.reactions) > 0:
						list(p.reactions)[0].delete(remove_orphans = True)
						for data in self.process_data.query(p.id):
							self.process_data.remove(data.id)

		for p in tqdm.tqdm(self.metabolites.query(re.compile('^protein_')), 'Pruning unnecessary ProcessedProtein reactions...', bar_format = bar_format):
			if isinstance(p, coralme.core.component.ProcessedProtein) and p.id not in skip:
				delete = True
				for rxn in p.reactions:
					if rxn.metabolites[p] < 0:
						delete = False
						break
				if delete:
					for rxn in list(p.reactions):
						self.process_data.remove(rxn.posttranslation_data.id)
						rxn.delete(remove_orphans = True)

		for p in tqdm.tqdm(self.metabolites.query(re.compile('^protein_')), 'Pruning unnecessary TranslatedGene reactions...', bar_format = bar_format):
			if isinstance(p, coralme.core.component.TranslatedGene) and p.id not in skip:
				delete = True
				for rxn in p.reactions:
					if rxn.metabolites[p] < 0 and not rxn.id.startswith('degradation'):
						delete = False
						break
				if delete:
					for rxn in p.reactions:
						p_id = p.id.replace('protein_', '')
						data = self.process_data.get_by_id(p_id)
						self.process_data.remove(data.id)
						rxn.delete(remove_orphans = True)

		removed_rna = set()
		for m in tqdm.tqdm(list(self.metabolites.query(re.compile('^RNA_'))), 'Pruning unnecessary TranscribedGene reactions...', bar_format = bar_format):
			delete = False if m.id in skip else True
			for rxn in m.reactions:
				if rxn.metabolites[m] < 0 and not rxn.id.startswith('DM_'):
					delete = False
			if delete:
				try:
					self.reactions.get_by_id('DM_' + m.id).remove_from_model(remove_orphans = True)
					if m in self.metabolites:
						# Defaults to subtractive when removing reaction
						m.remove_from_model()
				except KeyError:
					pass
				else:
					removed_rna.add(m.id)

		for t in tqdm.tqdm(self.reactions.query('transcription_TU'), 'Pruning unnecessary Transcriptional Units...', bar_format = bar_format):
			if t.id in skip:
				delete = False
			else:
				delete = True

			for product in t.products:
				if isinstance(product, coralme.core.component.TranscribedGene):
					delete = False

			t_process_id = t.id.replace('transcription_', '')
			if delete:
				t.remove_from_model(remove_orphans = True)
				self.process_data.remove(t_process_id)
			else:
				# gets rid of the removed RNA from the products
				self.process_data.get_by_id(t_process_id).RNA_products.difference_update(removed_rna)

			# update the TranscriptionReaction mRNA biomass stoichiometry with new RNA_products
			if not delete:
				t.update()

		return None

	def remove_genes_from_model(self, gene_list):
		for gene in tqdm.tqdm(gene_list, 'Removing gene(s) from ME-Model...', bar_format = bar_format):
			# defaults to subtractive when removing model
			self.metabolites.get_by_id('RNA_' + gene).remove_from_model()
			protein = self.metabolites.get_by_id('protein_'+gene)
			for cplx in protein.complexes:
				print('Complex \'{:s}\' removed from ME-Model.'.format(cplx.id))
				for rxn in cplx.metabolic_reactions:
					try:
						self.process_data.remove(rxn.id.split('_')[0])
					except ValueError:
						pass
					rxn.remove_from_model()

			protein.remove_from_model(destructive = True)

		# Remove all transcription reactions that now do not form a used transcript
		for tu in tqdm.tqdm(self.reactions.query('transcription_TU'), 'Removing unnecessary Transcriptional Units...', bar_format = bar_format):
			delete = True
			for product in tu.products:
				if isinstance(product, coralme.core.component.TranscribedGene):
					delete = False
			if delete:
				tu.remove_from_model(remove_orphans = True)
				t_process_id = tu.id.replace('transcription_', '')
				self.process_data.remove(t_process_id)

		return None

	def set_sasa_keffs(self, median_keff):
		# Get median SASA value considering all complexes in model
		sasa_list = []
		for met in tqdm.tqdm(self.metabolites, 'Processing Complexes...', bar_format = bar_format):
			cplx_sasa = 0.
			if not isinstance(met, coralme.core.component.Complex):
				continue
			cplx_sasa += met.formula_weight ** (3. / 4)
			sasa_list.append(cplx_sasa)
		median_sasa = numpy.median(numpy.array(sasa_list))

		# redo scaling average SASA to 65.
		for rxn in tqdm.tqdm(self.reactions, 'Processing Reactions...', bar_format = bar_format):
			if hasattr(rxn, 'keff') and rxn.complex_data is not None:
				sasa = rxn.complex_data.complex.formula_weight ** (3. / 4.)
				if sasa == 0:
					raise UserWarning('No SASA for reaction \'{:s}\'.'.format(rxn.id))
				rxn.keff = sasa * median_keff / median_sasa

		for data in tqdm.tqdm(self.process_data, 'Processing ProcessData...', bar_format = bar_format):
			sasa = 0.
			if isinstance(data, coralme.core.processdata.TranslocationData):
				continue
			if hasattr(data, 'keff') and hasattr(data, 'formula_weight') and data.enzyme is not None:
				cplxs = [data.enzyme] if type(data.enzyme) == str else data.enzyme
				for cplx in cplxs:
					sasa += self.metabolites.get_by_id(cplx).formula_weight ** (3. / 4)
				if sasa == 0:
					raise UserWarning('No SASA for reaction \'{:s}\'.'.format(data.id))
				data.keff = sasa * median_keff / median_sasa

		self.update()

		return None

	def update(self):
		new = []
		for r in self.reactions:
			if hasattr(r, 'update'):
				new.append(r)
		for r in tqdm.tqdm(new, 'Updating ME-Model Reactions...', bar_format = bar_format):
			_update(r)
		return None

	# me.update() cannot be paralelized without considering new constraints being added into the model.
	# New constraints must have a different name, so me.update() fails if two reactions are changed to add the same constraint:
	# ContainerAlreadyContains: Container '<optlang.container.Container object at 0x...>' already contains an object with name 'Name'.
	def _parallel_update(self):
		return None

	def get(self, x: typing.Union[cobra.core.object.Object, str]) -> cobra.core.object.Object:
		"""
		Return the element with a matching id from model.reactions or model.metabolites attributes.
		"""
		if isinstance(x, str) or isinstance(x, cobra.core.object.Object):
			if self.metabolites.has_id(x):
				return self.metabolites.get_by_id(x)
			elif self.reactions.has_id(x):
				return self.reactions.get_by_id(x)
			else:
				return
		else:
			return NotImplemented

	def query(self, x):
		"""
		Return the elements with a matching substring from model.reactions, model.metabolites, and model.process_data attributes.
		"""
		res = []
		x = x.replace('(', '\(').replace(')', '\)')
		res.append(self.metabolites.query(x))
		res.append(self.reactions.query(x))
		res.append(self.process_data.query(x))
		res = [ x for y in res for x in y ]
		if len(res) == 1:
			return res[0]
		else:
			return res

	# Originally developed by JDTB@UCSD, 2022
	def _optimize(self,
		max_mu = 1., min_mu = 0., fixed_mu: float = False,
		precision = 1e-6, solver_precision = 'quad', maxIter = 100,
		verbosity = False):
		"""
		Optimize the ME-Model using flux balance analysis.
		"""

		me = self
		# "repair"
		for met in me.metabolites:
			met._constraint_sense = 'E'

		if fixed_mu:
			from qminospy.me2 import ME_NLP
			me_nlp = ME_NLP(me)
			x, status, hs = me_nlp.solvelp(fixed_mu)
			me.solution.status = status
			me.solution.x_dict = { r:f for r,f in zip(me.reactions, x) }

		else:
			from qminospy.me1 import ME_NLP1
			# The object containing solveME methods--composite that uses a ME model object
			me_nlp = ME_NLP1(me, growth_key = self.mu)
			# Use bisection for now (until the NLP formulation is worked out)
			muopt, hs, xopt, cache = me_nlp.bisectmu(
				mumax = max_mu, mumin = min_mu,
				precision = precision, solver_precision = solver_precision, maxIter = maxIter,
				verbosity = verbosity)

		if me.solution:
			return True
		else:
			return False

	def optimize(self,
		max_mu = 1., min_mu = 0., fixed_mu: float = False,
		precision = 1e-6, solver_precision = 'quad', maxIter = 100,
		verbosity = False):

		me = self

		from coralme.solver.solver import ME_NLP
		me_nlp = ME_NLP(me)
		muopt, hs, xopt, cache = me_nlp.bisectmu(
				mumax = max_mu, mumin = min_mu,
				precision = precision, solver_precision = solver_precision, maxIter = maxIter,
				verbosity = verbosity)

		if me.solution:
			return True
		else:
			return False

	def relax_bounds(self):
		for rxn in self.reactions:
			if rxn.id == 'biomass_dilution':
				continue
			if hasattr(rxn.upper_bound, 'subs') or rxn.upper_bound > 0:
				rxn.upper_bound = 1000
			else:
				rxn.upper_bound = 0

			if hasattr(rxn.lower_bound, 'subs') or rxn.lower_bound > 0: # Is this OK?
				rxn.lower_bound = 0
			elif rxn.lower_bound < 0:
				rxn.lower_bound = -1000
