#!/usr/bin/python3
import copy
import pandas
import coralme
import multiprocessing
import itertools

# parallelization
from cobra.util import ProcessPool

def close_sink_and_solve(rxn_id):
	global _model
	global _muf
	model = copy.deepcopy(_model)
	model.reactions.get_by_id(rxn_id).bounds = (0, 0)
	model.optimize(max_mu = _muf, maxIter = 1)
	if not model.solution:
		return (rxn_id, False)
	else:
		return (rxn_id, True)

# Inspired from COBRApy parallelization
def _init_worker(me, muf):
    global _model
    global _muf
    _model = me
    _muf = muf

# Originally developed by JDTB@UCSD, 2022
def process_model(model, growth_key = 'mu'):
	dct = {}
	for met in model.metabolites:
		t = { 'c' : set(), 'p' : set() }
		#seen = [] #?
		for rxn in met.reactions:
			if rxn.id.startswith('BIOMASS_'):
				continue

			lb, ub = rxn.lower_bound, rxn.upper_bound

			# Replace 'growth_key' if model is a ME-Model
			if hasattr(lb, 'subs'):
				lb = lb.subs(growth_key, 1.)
			if hasattr(ub, 'subs'):
				ub = ub.subs(growth_key, 1.)

			coeff = rxn.metabolites[met]
			if hasattr(coeff, 'subs'):
				coeff = coeff.subs(growth_key, 1.)

			pos = 1 if coeff > 0 else -1
			rev = 1 if lb < 0 else 0
			fwd = 1 if ub > 0 else 0
			if pos*fwd == -1 or pos*rev == +1:
				t['c'].add(rxn.id)
			if pos*fwd == +1 or pos*rev == -1:
				t['p'].add(rxn.id)
		dct[met.id] = t
	return dct

def find_gaps(model, growth_key = 'mu'):
	g = {}
	dct = process_model(model, growth_key = growth_key)
	for met, t in dct.items():
		# not producing, not consuming, not uerever
		g[met] = { 'p' : 0, 'c' : 0, 'u' : 0 }
		if not t['c']:
			g[met]['c'] = 1
		if not t['p']:
			g[met]['p'] = 1
		if len(t['c']) == 1 and t['c'] == t['p']:
			g[met]['u'] = 1
	df = pandas.DataFrame.from_dict(g).T
	df = df[df.any(axis = 1)]
	df = df.sort_index()
	return df

def add_exchange_reactions(me, metabolites, prefix = 'SK_'):
	for met in metabolites:
		rxn_id = prefix + met
		if rxn_id not in me.reactions:
			r = coralme.core.reaction.MEReaction(rxn_id)
			me.add_reaction(r)
			r.add_metabolites({ met: -1 })
		else:
			r = me.reactions.get_by_id(rxn_id)
		r.bounds = (-1000, 1000)
		#print(r.id,r.lower_bound,r.upper_bound,r.reaction)
	return me

def brute_force_check(me, metabolites_to_add, objective_function = 'biomass_dilution', muf = 0.1):
#	 me.objective = objective_function
	print('  '*6 + 'Adding sink reactions for {:d} metabolites'.format(len(metabolites_to_add)))
	add_exchange_reactions(me, metabolites_to_add)
	print('  '*6 + 'Objective reaction: {:s}'.format(objective_function))

	me.optimize(max_mu = muf, maxIter = 1.)
	if not me.solution:
		print('  '*6 + 'No production capacity of objective')
		return 1

	print('  '*6 + 'Initial objective function value of {:f}'.format(me.solution.objective_value))
	rxns = []
	for idx, flux in me.solution.fluxes.items():
		if idx.startswith('SK_') and idx.split('SK_')[1] in metabolites_to_add:
			if abs(flux) > 0:
				rxns.append(idx)
			else:
				#print('Closing {}'.format(idx))
				me.reactions.get_by_id(idx).bounds = (0, 0)

	print('  '*6 + 'Sink reactions shortlisted to {:d} metabolites:'.format(len(rxns)))

	gaps = []
	#print()
	#print('  '*6 + 'Tested reaction Feasible? Progress Gaps found')
	#print('  '*6 + '=============== ========= ======== ==========')
	#for idx, rxn_id in enumerate(sorted(rxns)):
		#ex_rxn = me.reactions.get_by_id(rxn_id)
		#lb, ub = ex_rxn.bounds # save original bounds for later
		#ex_rxn.bounds = (0, 0)
		#me.optimize(max_mu = muf, maxIter = 1)
		#if not me.solution:
			#ex_rxn.bounds = (lb, ub)
			#gaps.append(rxn_id)
			#feasible = False
			#obj = 'gap'
		#else:
			#feasible = True
			#obj = me.solution.objective_value

		#print('  '*6 + '{:s}\t{:s}\t{:d}/{:d}\t{:d}'.format(rxn_id, str(feasible), idx+1, len(rxns), len(gaps)))

	with ProcessPool(processes = multiprocessing.cpu_count() - 2, initializer = _init_worker, initargs = (me, muf, )) as pool:
		for rxn_id, feasible in pool.imap_unordered(close_sink_and_solve, rxns):
			gaps.append((rxn_id, feasible))

	#return [ i.split('SK_')[1] for i in gaps ]
	return [ i.split('SK_')[1] for i,j in gaps if j ]
