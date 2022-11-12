#!/usr/bin/python3
import pickle
import logging

import coralme
from .helper_functions import find_gaps, brute_force_check, add_exchange_reactions

import multiprocessing

class METroubleshooter(object):
	"""METroubleshooter class for troubleshooting growth in a ME-Model

	This class contains methods to identify gaps to allow growth in
	a ME-Model. ME-Model must be saved as me_model.pickle in the folder
	containing the main organism.

	Parameters
	----------
	ME-Model : coralme.core.model.MEModel
		MEBuilder object containing the final version of the class
		Organism.


	"""
	def __init__(self, me_model):
		self.me_model = me_model
		self.me_model.curation_notes = { 'troubleshoot' : [] }

	def gap_find(self):
		#from draft_coralme.util.helper_functions import find_gaps
		sep = '~ '*6
		print(sep+'Finding gaps from the M-Model only...')
		m_gaps = find_gaps(self.me_model.gem)

		print(sep+'Finding gaps in the ME-Model...')
		me_gaps = find_gaps(self.me_model, growth_key = self.me_model.mu)

		idx = list(set(me_gaps.index) - set(m_gaps.index))
		new_gaps = me_gaps.loc[idx]

		filt1 = new_gaps['p'] == 1
		filt2 = new_gaps['c'] == 1
		filt3 = ~(new_gaps['u'] == 1)

		deadends = list(new_gaps[filt1 | filt2 | filt3].index)
		deadends = [ x for x in deadends if 'biomass' not in x ]
		deadends = sorted(deadends)
		print('  '*6 + '{:d} metabolites were identified as deadends.'.format(len(deadends)))
		for met in deadends:
			name = self.me_model.metabolites.get_by_id(met).name
			print('  '*6 + '{:s}: {:s}'.format(met, 'No name' if name == '' else name))
		return deadends

	def gap_fill(self, deadends = [], met_type = 'Metabolite', mu_test = 0.1):
		sep = '~ '*6
		#from draft_coralme.util.helper_functions import add_exchange_reactions
		if deadends:
			print(sep+'Adding a sink reaction for each identified deadend metabolite.')
			add_exchange_reactions(self.me_model, deadends)

		print(sep+'Optimizing ME-Model...')
		self.me_model.optimize(max_mu = 0.1, maxIter = 1)

		bf_gaps = []
		works = False

		if self.me_model.solution:
			print('  '*6 + 'Gapfilled ME-Model is feasible with growth rate {:f}.'.format(self.me_model.solution.objective_value))
			works = True

		else:
			#from draft_coralme.util.helper_functions import brute_force_check
			print('  '*6 + 'Provided set of sink reactions for deadend metabolites does not allow growth.')
			print(sep+'Step 3. Looking for gaps in the network using brute force.')

			if isinstance(met_type, str):
				mets = [ met.id for met in self.me_model.metabolites \
					if isinstance(met, getattr(coralme.core.component, met_type)) and not met.id.endswith('_e') ]
			elif isinstance(met_type, list):
				mets = met_type
			else:
				raise ValueError('Argument of type \'{:s}\' is not valid. It must be a string or a list of strings.'.format(type(met_type)))

			# remove metabolites that are feed into the model through transport reactions
			medium = set([ '{:s}_c'.format(x[3:-2]) for x in self.me_model.gem.medium.keys() ])
			mets = set(mets).difference(medium)

			# filter manually
			mets = set(mets).difference(set(['ppi_c', 'ACP_c']))
			mets = set(mets).difference(set(['ade_c']))
			mets = set(mets).difference(set(['adp_c', 'amp_c', 'atp_c']))
			mets = set(mets).difference(set(['cdp_c', 'cmp_c', 'ctp_c']))
			mets = set(mets).difference(set(['gdp_c', 'gmp_c', 'gtp_c']))
			mets = set(mets).difference(set(['udp_c', 'ump_c', 'utp_c']))
			mets = set(mets).difference(set(['dadp_c', 'dcdp_c', 'dgdp_c', 'dtdp_c', 'dudp_c']))
			mets = set(mets).difference(set(['damp_c', 'dcmp_c', 'dgmp_c', 'dtmp_c', 'dump_c']))
			mets = set(mets).difference(set(['datp_c', 'dctp_c', 'dgtp_c', 'dttp_c', 'dutp_c']))
			mets = set(mets).difference(set(['nad_c', 'nadh_c', 'nadp_c', 'nadph_c']))
			mets = set(mets).difference(set(['5fthf_c', '10fthf_c', '5mthf_c', 'dhf_c', 'methf_c', 'mlthf_c', 'thf_c']))
			mets = set(mets).difference(set(['fad_c', 'fadh2_c', 'fmn_c']))
			mets = set(mets).difference(set(['coa_c']))

			bf_gaps = brute_force_check(self.me_model, mets, muf = mu_test)
			if bf_gaps == 1:
				works = False
			else:
				works = True

		return bf_gaps, works

	def troubleshoot(self, mu_test = 0.001):
		sep = '~ '*6

		print(sep+'Troubleshooting started...')
		print(sep+'Step 1. Find topological gaps in the ME-Model...')
		deadends = self.gap_find()

		if len(deadends) != 0:
			print(sep+'Step 2. Solve gap-filled ME-Model with all identified deadend metabolites...')

		works = False
		bf_gaps = 0
		mods = []

		self.me_model.optimize(max_mu = mu_test, maxIter = 1)

		if self.me_model.solution:
			print('  '*6 + 'Original ME-Model is feasible with test growth rate of {:f} 1/h'.format(mu_test))
			#works = True
			mods.append('clean')

		else:
			bf_gaps, works = self.gap_fill(deadends = deadends)

			if works:
				with open('MEModel-step3-{:s}-gapfilled.pkl'.format(self.me_model.id), 'wb') as outfile:
					pickle.dump(self.me_model, outfile)
				mods.append('gapfilled')

			else:
				print('  '*6 + 'Brute forcing with all metabolic metabolites did not allow growth.')
				print(sep+'Step 4. Trying different groups of E-matrix components')

				met_types = [ 'Complex', 'GenerictRNA', 'TranscribedGene', 'TranslatedGene', 'ProcessedProtein', 'GenericComponent' ]
				for mt in met_types:
					print('  '*6 + 'Trying type \'{:s}\' of metabolites...'.format(mt))
					#self.me_model = me_builder.load_me()
					self.me_model.relax_bounds()

					for rxn in self.me_model.reactions.query('biomass_to_biomass'):
						rxn.upper_bound = mu_test

					self.me_model.reactions.protein_biomass_to_biomass.lower_bound = mu_test/100 # Needed to enforce protein production
					deadends = [ met.id for met in self.me_model.metabolites if isinstance(met, coralme.core.component.Metabolite) ]

					bf_gaps, works = self.gap_fill(deadends = deadends, met_type = mt)
					if works:
						mods.append(mt)
						break

		if deadends:
			self.me_model.curation_notes['troubleshoot'].append({
				'msg':'Some deadends were identified',
				'triggered_by':deadends,
				'importance':'high',
				'to_do':'Fix the deadends by adding reactions or solving other warnings.'})

		if works: # Meaning it could grow in any case
			if isinstance(bf_gaps, list) and bf_gaps:
				self.me_model.curation_notes['troubleshoot'].append({
					'msg':'Some metabolites are necessary for growth',
					'triggered_by':bf_gaps,
					'importance':'critical',
					'to_do':'Fix the gaps by adding reactions or solving other warnings. If some items are from the E-matrix, fix these first!'})

			print(sep+'Final step. Fully optimizing with precision 1e-6 and saving ME-Model...')
			self.me_model.optimize(max_mu = 3.0, precision = 1e-6)
			print('  '*6 + 'Gapfilled ME-Model is feasible with growth rate {:f}.'.format(self.me_model.solution.objective_value))

			# delete added sink reactions with lb == 0 and up == 0
			for rxn in self.me_model.reactions.query('^SK_'):
				if rxn.lower_bound == 0 and rxn.upper_bound == 0:
					self.me_model.remove_reactions([rxn])

			with open('MEModel-step3-{:s}-{:s}.pkl'.format(self.me_model.id, '_'.join(mods)), 'wb') as outfile:
			#me.solution.fluxes.to_csv('{}_converged_fluxes.csv'.format('_'.join(mods)))
			#me_builder.save_me(me, 'me_model_{}_converged.pickle'.format('_'.join(mods)))
				pickle.dump(self.me_model, outfile)
		else:
			logging.warning('Troubleshooter failed to determine a set of problematic metabolites.')

		return None
