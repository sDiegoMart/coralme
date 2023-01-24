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
import coralme

def _update(MEReaction):
	"""updates all component reactions"""
	MEReaction.update()
	return None

class MEModel(cobra.core.model.Model):
	def __init__(self, name, mu = 'mu'):
		cobra.Model.__init__(self, name)

		self.global_info = {
			'domain' : 'Prokaryote',

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
			'include_pseudo_genes' : False,

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
				'c' : 'Cytoplasm',
				'e' : 'Extracellular',
				'p' : 'Periplasm',
				'mc': 'ME-model Constraint'
				},

			'START_tRNA' : [],
			'rna_components' : [],
			'knockouts' : [],
			'genome_mods' : {},
			'trna_misacylation' : {},
			'trna_to_codon' : {},

			'gam' : 45.,
			'ngam' : 1.,
			'unmodeled_protein_fraction' : 0.36,

			'braun\'s_lipoprotein' : [],
			'braun\'s_lipid_mod' : 'murein5px4p_p',
			'braun\'s_lpp_flux' : -0.0,
			'braun\'s_murein_flux' : -0.0,
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
		# Solved: Check if variable type is sympy.core.symbol.Symbol or float
		# Solved: Removed _populate_solver from reactions -> no need to modify optlang
		self._biomass_dilution.upper_bound = self.mu
		self._biomass_dilution.lower_bound = self.mu

		# Maintenance energy
		self._gam = self.global_info['gam'] # default value
		self._ngam = self.global_info['ngam'] # default value

		"""
		Unmodeled protein is handled by converting protein_biomass to
		biomass, and requiring production of the appropriate amount of dummy
		protein
		"""
		self._unmodeled_protein_fraction = self.global_info['unmodeled_protein_fraction'] # default value

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

	def remove_reactions(self, reactions, remove_orphans=False):
		"""Remove reactions from the model.

		Parameters
		----------
		reactions : list
			A list with reactions (`cobra.Reaction`), or their id's, to remove

		remove_orphans : bool
			Remove orphaned genes and metabolites from the model as well

		"""
		if isinstance(reactions, str) or hasattr(reactions, "id"):
			reactions = [reactions]

		for reaction in reactions:
			# Make sure the reaction is in the model
			try:
				reaction = self.reactions[self.reactions.index(reaction)]
			except ValueError:
				logging.warning("%s not in %s" % (reaction, self))

			else:
				#forward = reaction.forward_variable
				#reverse = reaction.reverse_variable

				#self.remove_cons_vars([forward, reverse])
				self.reactions.remove(reaction)
				reaction._model = None

				for met in reaction._metabolites:
					if reaction in met._reaction:
						met._reaction.remove(reaction)
						if remove_orphans and len(met._reaction) == 0:
							self.remove_metabolites(met)

				for gene in reaction._genes:
					if reaction in gene._reaction:
						gene._reaction.remove(reaction)
						if remove_orphans and len(gene._reaction) == 0:
							self.genes.remove(gene)

				# remove reference to the reaction in all groups
				associated_groups = self.get_associated_groups(reaction)
				for group in associated_groups:
					group.remove_members(reaction)

	def add_biomass_constraints_to_model(self, biomass_types):
		for biomass_type in tqdm.tqdm(biomass_types, 'Adding biomass constraint(s) into the ME-model...', bar_format = bar_format):
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
			logging.warning('Adding GAM (ATP requirement for growth) reaction into the ME-model.')
			self.add_reactions([coralme.core.reaction.SummaryVariable('GAM')])
			self.reactions.GAM.lower_bound = self.mu
			self.reactions.GAM.upper_bound = self.mu
		#atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1, 'pi_c': 1} # charges: -4, 0 => -3, +1, -2
		atp_hydrolysis = self.process_data.get_by_id('atp_hydrolysis').stoichiometry
		for met, coeff in atp_hydrolysis.items():
			self.reactions.GAM.add_metabolites({met: value * coeff}, combine = False)
		self._gam = value

		# check stoichiometry
		if self.reactions.GAM.check_mass_balance() == {'charge': -1.0, 'H': -1.0}:
			self.reactions.GAM._metabolites.update({self.metabolites.h_c : +1})

	@property
	def ngam(self):
		return self._ngam

	@ngam.setter
	def ngam(self, value):
		if 'ATPM' not in self.reactions:
			logging.warning('Adding ATPM (ATP requirement for maintenance) reaction into the ME-model.')
			#atp_hydrolysis = {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'h_c': 1, 'pi_c': 1} # charges: -4, 0 => -3, +1, -2
			atp_hydrolysis = self.process_data.get_by_id('atp_hydrolysis').stoichiometry
			self.add_reactions([coralme.core.reaction.SummaryVariable('ATPM')])
			self.reactions.ATPM.add_metabolites(atp_hydrolysis)
		self.reactions.ATPM.lower_bound = value
		self.reactions.ATPM.upper_bound = value
		self._ngam = value

		# check stoichiometry
		if self.reactions.ATPM.check_mass_balance() == {'charge': -1.0, 'H': -1.0}:
			self.reactions.ATPM._metabolites.update({self.metabolites.h_c : +1})

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

	# ME-model methods
	def get_metabolic_flux(self, solution = None):
		"""Extract the flux state for Metabolic reactions."""
		if solution is None:
			solution = self.solution
		if solution.status != 'optimal':
			raise ValueError('Solution status \'{:s}\' is not \'optimal\'.'.format(solution.status))
		flux_dict = {r.id: 0 for r in tqdm.tqdm(list(self.stoichiometric_data), 'Building reaction dictionary...', bar_format = bar_format)}
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-model Reactions...', bar_format = bar_format):
			if isinstance(reaction, coralme.core.reaction.MetabolicReaction):
				m_reaction_id = reaction.stoichiometric_data.id
				if reaction.reverse:
					flux_dict[m_reaction_id] -= solution.fluxes[reaction.id]
				else:
					flux_dict[m_reaction_id] += solution.fluxes[reaction.id]
			# SummaryVariable in M-model
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
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-model Reactions...', bar_format = bar_format):
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
		for reaction in tqdm.tqdm(self.reactions, 'Processing ME-model Reactions...', bar_format = bar_format):
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
			for met, value in rxn.metabolites.items():
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
		for gene in tqdm.tqdm(gene_list, 'Removing gene(s) from ME-model...', bar_format = bar_format):
			# defaults to subtractive when removing model
			self.metabolites.get_by_id('RNA_' + gene).remove_from_model()
			protein = self.metabolites.get_by_id('protein_'+gene)
			for cplx in protein.complexes:
				print('Complex \'{:s}\' removed from ME-model.'.format(cplx.id))
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
		for r in tqdm.tqdm(new, 'Updating ME-model Reactions...', bar_format = bar_format):
			_update(r)
		return None

	# me.update() cannot be paralelized without considering new constraints being added into the model.
	# New constraints must have a different name, so me.update() fails if two reactions are changed to add the same constraint:
	# ContainerAlreadyContains: Container '<optlang.container.Container object at 0x...>' already contains an object with name 'Name'.
	def _parallel_update(self):
		return NotImplemented

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

	def construct_lp_problem(self, keys = False):
		# populate empty dictionaries with stoichiometry
		Sf = dict()
		Se = dict()

		#types = (sympy.core.add.Add, sympy.core.mul.Mul)

		# check how many variables are in the ME-model
		atoms = set()

		for idx, rxn in enumerate(self.reactions):
			for met, value in rxn.metabolites.items():
				met_index = self.metabolites.index(met)
				if hasattr(value, 'subs'):
					atoms.add(list(value.free_symbols)[0])
					Se[met_index, idx] = value
				else:
					Sf[met_index, idx] = value

		lb, ub = zip(*[ rxn.bounds for rxn in self.reactions ])
		b = [ m._bound for m in self.metabolites ] # accumulation
		#c = [ r.objective_coefficient for r in self.reactions ]
		c = [ 0 for r in self.reactions ]
		# constraint sense eventually be in the metabolite...
		cs = [ 'E' for m in self.metabolites ]

		if keys:
			# replace symbols before calling the solver
			Se = { k:float(x.xreplace(keys)) if hasattr(x, 'subs') else x for k,x in Se.items() }
			lb = [ float(x.xreplace(keys)) if hasattr(x, 'subs') else x for x in lb ]
			ub = [ float(x.xreplace(keys)) if hasattr(x, 'subs') else x for x in ub ]

			Sf.update(Se)

		return Sf, Se, lb, ub, b, c, cs, atoms

	def optimize(self,
		max_mu = 1., min_mu = 0., maxIter = 100, lambdify = True,
		tolerance = 1e-6, precision = 'quad', verbose = False):

		# check options
		tolerance = tolerance if tolerance >= 1e-15 else 1e-6
		precision = precision if precision in [ 'quad', 'double', 'dq', 'dqq' ] else 'quad'

		# populate with stoichiometry, no replacement of mu's
		Sf, Se, lb, ub, b, c, cs, atoms = self.construct_lp_problem()

		if len(atoms) > 1:
			print('Use `me_model.map_feasibility()` to obtain the boundary of feasibility solutions.')
			print('Optimization will proceed replacing all growth keys with the same value.')

		if lambdify:
			fn = numpy.vectorize(lambda x: sympy.lambdify(list(atoms), x))
			lambdas = { k:v for k,v in zip(Se.keys(), fn(list(Se.values()))) }
		else:
			lambdas = None

		from coralme.solver.solver import ME_NLP
		#me_nlp = ME_NLP(me)
		me_nlp = ME_NLP(Sf, Se, b, c, lb, ub, cs, atoms, lambdas)
		muopt, xopt, yopt, basis, stat = me_nlp.bisectmu(
				mumax = max_mu,
				mumin = min_mu,
				maxIter = maxIter,
				tolerance = tolerance,
				precision = precision,
				verbose = verbose)

		#f = sum([ rxn.objective_coefficient * xopt[idx] for idx, rxn in enumerate(self.reactions) ])
		x_primal = xopt[ 0:len(self.reactions) ]   # The remainder are the slacks
		x_dict = { rxn.id : xopt[idx] for idx, rxn in enumerate(self.reactions) }
		#y = pi
		# J = [S; c]
		y_dict = { met.id : yopt[idx] for idx, met in enumerate(self.metabolites) }
		#y_dict['linear_objective'] = y[len(y)-1]

		if stat == 'optimal':
			#self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)
			self.solution = cobra.core.Solution(
				objective_value = muopt,
				status = stat,
				fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
				reduced_costs = None,
				shadow_prices = None,
				)
			return True
		else:
			return False

	def feasibility(self, keys = { sympy.Symbol('mu', positive = True) : 1. }, tolerance = 1e-6, precision = 'quad', **kwargs):
		# check options
		tolerance = tolerance if tolerance >= 1e-15 else 1e-6
		precision = precision if precision in [ 'quad', 'double', 'dq', 'dqq' ] else 'quad'
		for key in list(keys.keys()):
			if isinstance(key, sympy.Symbol):
				pass
			else:
				keys[sympy.Symbol(key, positive = True)] = keys.pop(key)

		# populate with stoichiometry with replacement of mu's
		Sf, Se, lb, ub, b, c, cs, atoms = kwargs.get('lp', self.construct_lp_problem(keys = keys))

		from coralme.solver.solver import ME_NLP
		#me_nlp = ME_NLP(me)
		me_nlp = ME_NLP(Sf, dict(), b, c, lb, ub, cs, set(keys.keys()), None)
		muopt, xopt, yopt, basis, stat = me_nlp.bisectmu(
				mumax = 1., # mu was already replaced and maxIter is one, so a value here doesn't matter
				mumin = 0.,
				maxIter = 1,
				tolerance = tolerance,
				precision = precision,
				verbose = False)

		#f = sum([ rxn.objective_coefficient * xopt[idx] for idx, rxn in enumerate(self.reactions) ])
		x_primal = xopt[ 0:len(self.reactions) ]   # The remainder are the slacks
		x_dict = { rxn.id : xopt[idx] for idx, rxn in enumerate(self.reactions) }
		#y = pi
		# J = [S; c]
		y_dict = { met.id : yopt[idx] for idx, met in enumerate(self.metabolites) }
		#y_dict['linear_objective'] = y[len(y)-1]

		if stat == 'optimal':
			#self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)
			self.solution = cobra.core.Solution(
				objective_value = list(keys.values())[0],
				status = stat,
				fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
				reduced_costs = None,
				shadow_prices = None,
				)
			return True
		else:
			return False

	def map_feasibility(self, keys = { sympy.Symbol('mu', positive = True) : 1. }, tolerance = 1e-6, precision = 'quad'):
		return NotImplemented

	# Originally developed by JDTB@UCSD, 2022
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
