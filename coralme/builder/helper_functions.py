#!/usr/bin/python3
import os
import re
import sys
import copy
import json

import logging
log = logging.getLogger(__name__)

import sympy
import pandas
import tqdm

from ast import parse as ast_parse, Name, And, Or, BitOr, BitAnd, BoolOp, Expression, NodeTransformer

import cobra
import coralme

import collections

# from cobrame without changes
def get_base_complex_data(model, complex_id):
	"""
	If a complex is modified in a metabolic reaction it will not
	have a formation reaction associated with it. This function returns
	the complex data of the "base" complex, which will have the subunit
	stoichiometry of that complex
	"""

	# First try unmodified complex id
	try_1 = complex_id.split('_mod_')[0]
	if try_1 in model.process_data:
		return model.process_data.get_by_id(try_1)

	count = 0
	try_2 = complex_id.split('_mod_')[0] + '_'
	for i in model.process_data.query(try_2):
		if isinstance(i, coralme.core.processdata.ComplexData):
			count += 1
			data = i

	if count == 0:
		raise UserWarning('No base complex found for \'{:s}\'.'.format(complex_id))

	elif count > 1:
		raise UserWarning('More than one possible base complex found for \'{:s}\'.'.format(complex_id))

	return data

# from cobra with changes
def parse_composition(tmp_formula) -> dict:
	element_re = re.compile("([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?)")

	"""Break the chemical formula down by element."""
	#tmp_formula = self.formula
	# commonly occuring characters in incorrectly constructed formulas
	if "*" in tmp_formula:
		warn(f"invalid character '*' found in formula '{self.formula}'")
		tmp_formula = self.formula.replace("*", "")
	if "(" in tmp_formula or ")" in tmp_formula:
		warn(f"parenthesis found in formula '{self.formula}'")
		return
	composition = {}
	parsed = element_re.findall(tmp_formula)
	for element, count in parsed:
		if count == "":
			count = 1
		else:
			try:
				count = float(count)
				int_count = int(count)
				if count == int_count:
					count = int_count
				else:
					warn(f"{count} is not an integer (in formula {self.formula})")
			except ValueError:
				warn(f"failed to parse {count} (in formula {self.formula})")
				#self.elements = {}
				return
		if element in composition:
			composition[element] += count
		else:
			composition[element] = count
	return composition

# Originally developed by JDTB, UCSD (2022)
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

def change_reaction_id(model,old_id,new_id):
	old_rxn = model.reactions.get_by_id(old_id)
	rxn = cobra.Reaction(new_id)
	model.add_reactions([rxn])
	for k,v in old_rxn.metabolites.items():
		rxn.add_metabolites({k.id:v})
	rxn.bounds = old_rxn.bounds
	rxn.name = old_rxn.name
	rxn.subsystem = old_rxn.subsystem
	rxn.notes = old_rxn.notes
	rxn.gene_reaction_rule = old_rxn.gene_reaction_rule
	model.remove_reactions([old_rxn])

def listify_gpr(expr,level = 0):
	if level == 0:
		return listify_gpr(cobra.core.gene.parse_gpr(str(expr))[0], level = 1)
	if isinstance(expr, cobra.core.gene.GPR):
		return listify_gpr(expr.body, level = 1) if hasattr(expr, "body") else ""
	elif isinstance(expr, Name):
		return expr.id
	elif isinstance(expr, BoolOp):
		op = expr.op
		if isinstance(op, Or):
			str_exp = list([listify_gpr(i, level = 1) for i in expr.values])
		elif isinstance(op, And):
			str_exp = tuple([listify_gpr(i, level = 1) for i in expr.values])
		return str_exp
	elif expr is None:
		return ""
	else:
		raise TypeError("unsupported operation  " + repr(expr))

def get_combinations(l_gpr):
	l_gpr = [[i] if isinstance(i,str) else i for i in l_gpr]
	return list(product(*l_gpr))

def get_chain(l_gpr):
	l_gpr = [[i] if isinstance(i,str) else i for i in l_gpr]
	return list(chain(*l_gpr))

def _expand_gpr(l_gpr):
	if isinstance(l_gpr,str):
		return l_gpr
	elif isinstance(l_gpr,list):
		return ' or '.join([expand_gpr(i) for i in l_gpr])
	elif isinstance(l_gpr,tuple):
		if len(l_gpr) == 1:
			return expand_gpr(l_gpr[0])
		elif any(isinstance(i,list) for i in l_gpr):
			return ' or '.join(chain([expand_gpr(i) for i in get_combinations(l_gpr)]))
		else:
			return ' and '.join(get_chain(l_gpr))

# def generify_gpr(l_gpr,rxn_id,d={}):
# 	if isinstance(l_gpr,str):
# 		name = l_gpr
# 		return name,d
# 	elif isinstance(l_gpr,list):
# 		l = []
# 		for i in l_gpr:
# 			n,d = generify_gpr(i,rxn_id,d=d)
# 			l.append(n)
# 		base_name = 'generic_{}'.format(rxn_id)
# 		name = '{}_{}'.format(base_name,len([i for i in d if base_name in i]))
# 		d[name] = ' or '.join(l)
# 		return name,d
# 	elif isinstance(l_gpr,tuple):
# 		l = []
# 		for i in l_gpr:
# 			n,d = generify_gpr(i,rxn_id,d=d)
# 			l.append(n)
# 		base_name = 'CPLX_{}'.format(rxn_id)
# 		name = '{}-{}'.format(base_name,len([i for i in d if base_name in i]))
# 		d[name] = ' and '.join(l)
# 		return name,d

def print_check(i,l_gpr,T):
	print(i)
	print(l_gpr)
	print(T)
	print()

def get_tree(l_gpr,T={}):
	if isinstance(l_gpr,str):
		return l_gpr
	else:
		if isinstance(l_gpr,list):
			op = 'or'
		elif isinstance(l_gpr,tuple):
			op = 'and'
		T[op] = []
		for idx,i in enumerate(l_gpr):
			d = {}
			T[op].append(get_tree(i,T=d))
		return T

def append_graph(G,g):
	if G == '$':
		return g.copy()
	if isinstance(G,dict):
		for k,v in G.items():
			G[k] = append_graph(v,g)
		return G
def concatenate_graphs(L,r=[],threshold=100):
	if L == "STOP":
		return L
	if r:
		for i in r:
			if get_size(L) > threshold:
				return "STOP"
			L = append_graph(L,i)
		return L
	elif isinstance(L,list):
		if len(L) == 1:
			return L[0]
		else:
			b = L[0]
			r = L[1:]
			L = concatenate_graphs(b,r,threshold=threshold)
		return L

def get_size(G):
	return len(re.findall("\$",str(G)))

def get_graph(T,G={},threshold=100):
	if isinstance(T,str):
		if T in G:
			T = T + '_REPETITIONMARK_' + str(len(G))
		G[T] = '$'
		return G
	elif isinstance(T,dict):
		if 'and' in T:
			l = []
			for i in T['and']:
				d = {}
				g = get_graph(i,d,threshold=threshold)
				if g == "STOP":
					return g
				l.append(g)
			d = concatenate_graphs(l,threshold=threshold)
			if d == "STOP":
				return d
			for k,v in d.items():
				if k in G:
					k = k + '_REPETITIONMARK_' + str(len(G))
				G[k] = v
			return G
		elif 'or' in T:
			for i in T['or']:
				G = get_graph(i,G,threshold=threshold)
				if G == "STOP":
					return G
		return G

def traverse_graph(G,L = [], C = []):
	if G == '$':
		C.append(L)
		return L,C
	if isinstance(G,dict):
		for k,v in G.items():
			k = k.split('_REPETITIONMARK_')[0]
			l = L + [k]
			l,C = traverse_graph(v,l,C)
		return L,C

def expand_gpr(rule,threshold=100):
	l = listify_gpr(rule)
	T = get_tree(l,T={})
	G = get_graph(T,G={},threshold=threshold)
	if G == "STOP":
		return G
	return traverse_graph(G,L=[],C=[])[1]

def generify_gpr(l_gpr,rxn_id,d={},generic_gene_dict={}):
	if isinstance(l_gpr,str):
		name = l_gpr
		return name,d
	elif isinstance(l_gpr,list):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d,generic_gene_dict=generic_gene_dict)
			l.append(n)
		existing_generic = find_match(generic_gene_dict,l)
		if existing_generic:
			name = existing_generic
		else:
			base_name = 'generic_{}'.format(rxn_id)
			name = '{}_{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' or '.join(l)
		return name,d
	elif isinstance(l_gpr,tuple):
		l = []
		for i in l_gpr:
			n,d = generify_gpr(i,rxn_id,d=d,generic_gene_dict=generic_gene_dict)
			l.append(n)
		base_name = 'CPLX_{}'.format(rxn_id)
		name = '{}-{}'.format(base_name,len([i for i in d if base_name in i]))
		d[name] = ' and '.join(l)
		return name,d

def process_rule_dict(n,rule_dict,gene_dict,protein_mod):
	corrected_ids = {}
	for cplx,rule in rule_dict.items():
		cplx_id = 0
		if 'CPLX' in cplx:
			rule_gene_list = rule.split(" and ")
			identified_genes = rule_gene_list
			cplx_id = find_match(gene_dict,identified_genes)
		if not cplx_id:
			cplx_id = cplx
		corrected_ids[cplx] = cplx_id
	corrected_rule_dict = {}

	for cplx,rule in rule_dict.items():
		if cplx in corrected_ids:
			cplx_id = corrected_ids[cplx]
		else:
			cplx_id = cplx
		# if cplx_id in protein_mod
		if cplx_id in protein_mod["Core_enzyme"].values:
			cplx_mod_id = protein_mod[
				protein_mod["Core_enzyme"].str.contains(cplx_id)
			].index[0]
			if "Oxidized" in cplx_mod_id:
				cplx_mod_id = cplx_mod_id.split("_mod_Oxidized")[0]
			if corrected_ids[n] == cplx_id:
				rule = corrected_ids.pop(n)
				corrected_ids[n] = cplx_mod_id
			cplx_id = cplx_mod_id
		for c,cid in corrected_ids.items():
			regex = '{}(?!\d)'
			corrected_rule_dict[cplx_id] = re.sub(regex.format(c), cid, rule)
			rule = corrected_rule_dict[cplx_id]
	return corrected_ids[n],corrected_rule_dict

def find_match(d,items):
    for c, cg in d.items():
        if not cg: continue
        if isinstance(cg,str):
            cg = [re.findall('.*(?=\(\d*\))', g)[0] for g in cg.split(' AND ')]
        if set(cg) == set(items):
            return c
    return 0

# Originally developed by JDTB@UCSD, 2022
# Modified by RSP@UCSD, 2022
def process_model(model, growth_key = sympy.Symbol('mu', positive = True)):
	dct = {}
	for met in model.metabolites:
		filter1 = type(met) == cobra.core.metabolite.Metabolite or type(met) == coralme.core.component.Metabolite
		filter2 = met.id.startswith('trna')
		filter3 = met.id.endswith('trna_c')

		if filter1 and not filter2 and not filter3:
			t = { 'c' : set(), 'p' : set() }
			#seen = [] #?
			for rxn in met.reactions:
				if rxn.id.startswith('BIOMASS_'):
					continue

				lb, ub = rxn.lower_bound, rxn.upper_bound

				# Replace 'growth_key' if model is a ME-model
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

def find_gaps(model, growth_key = sympy.Symbol('mu', positive = True)):
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
	rxns = []
	for met in metabolites:
		rxn_id = prefix + met
		if rxn_id not in me.reactions:
			r = coralme.core.reaction.MEReaction(rxn_id)
			me.add_reactions([r])
			r.add_metabolites({ met: -1 })
		else:
			r = me.reactions.get_by_id(rxn_id)
		r.bounds = (-10, 1000)
		rxns.append(r)
		#print(r.id,r.lower_bound,r.upper_bound,r.reaction)
	return rxns

def exchange_single_model(me, flux_dict = 0, solution=0):
	complete_dict = {'id':[],'name':[],'reaction':[],'lower_bound':[],'upper_bound':[],'flux':[]}

	if solution:
		flux_dict = solution.fluxes
	elif not flux_dict:
		flux_dict = me.solution.fluxes

	for rxn in me.reactions:
		try:
			if rxn.reactants and rxn.products:
				continue
		except:
			continue
		flux = flux_dict[rxn.id]

		if not flux:
			continue
		rxn_name = rxn.name
		reaction = rxn.reaction
		lb = rxn.lower_bound
		ub = rxn.upper_bound

		complete_dict['id'].append(rxn.id)
		complete_dict['name'].append(rxn_name)
		complete_dict['reaction'].append(reaction)
		complete_dict['lower_bound'].append(lb)
		complete_dict['upper_bound'].append(ub)
		complete_dict['flux'].append(flux)


	df = pandas.DataFrame(complete_dict).set_index('id')
	return df

def get_metabolites_from_pattern(model,pattern):
    met_list = []
    for met in model.metabolites:
        if pattern in met.id:
            met_list.append(met.id)
    return met_list

def get_met_coeff(stoich,growth_rate,growth_key='mu'):
	if hasattr(stoich, 'subs'):
		try:
			return float(stoich.subs(growth_key,growth_rate))
		except:
			return None
	return stoich

def summarize_reactions(model,met_id,only_types=(),ignore_types = ()):
	reactions = get_reactions_of_met(model,met_id,only_types=only_types,
								 ignore_types=ignore_types,verbose=False)
	d = {}
	for r in reactions:
		if r.bounds == (0,0):
			continue
		d[r.id] = {
			'name':r.name,
			'gene_reaction_rule':r.gene_reaction_rule,
			'reaction':r.reaction,
			'notes':r.notes if r.notes else ''
		}
	df = pandas.DataFrame.from_dict(d).T
	return df[['name','gene_reaction_rule','reaction','notes']] if not df.empty else 'No reaction found'


def flux_based_reactions(model,
						 met_id,
						 only_types=(),
						 ignore_types = (),
						 threshold = 0.,
						 flux_dict=0,
						 solution = None,
						 keffs=False,
						 verbose=False):

	if not flux_dict:
		#flux_dict = model.solution.x_dict
		if not hasattr(model,'solution') or not model.solution:
			if solution is not None:
				flux_dict = solution.fluxes
			else:
				print('No solution in model object')
				flux_dict = {r.id:0. for r in model.reactions}
		else:
			flux_dict = model.solution.fluxes
	mu = model.mu if hasattr(model,'mu') else ''
	reactions = get_reactions_of_met(model,met_id,only_types=only_types,
									 ignore_types=ignore_types,verbose=False,growth_key=mu)
	if len(reactions) == 0:
		print('No reactions found for {}'.format(met_id))
		return

	met = model.metabolites.get_by_id(met_id)
	result_dict = {}
	g = flux_dict.get('biomass_dilution',None)
	for rxn in (tqdm.tqdm(reactions) if verbose else reactions):
	#for rxn in tqdm.tqdm(reactions):
		f = flux_dict[rxn.id]
		result_dict[rxn.id] = {}
		if f:
			coeff = get_met_coeff(rxn.metabolites[met],
								  g,
								growth_key=model.mu if hasattr(model,"mu") else None)

		else:
			coeff = 0
		if coeff is None:
			print('Could not convert expression to float in {}'.format(rxn.id))
			continue
		try:
			result_dict[rxn.id]['lb'] = rxn.lower_bound if isinstance(rxn.lower_bound, sympy.Symbol) else float(rxn.lower_bound)
			result_dict[rxn.id]['ub'] = rxn.upper_bound if isinstance(rxn.upper_bound, sympy.Symbol) else float(rxn.upper_bound)
		except:
			print('Could not convert bounds to float in {}'.format(rxn.id))
		result_dict[rxn.id]['rxn_flux'] = f
		result_dict[rxn.id]['met_flux'] = f*coeff
		result_dict[rxn.id]['reaction'] = rxn.reaction
		if keffs:
			result_dict[rxn.id]['keff'] = rxn.keff if hasattr(rxn,'keff') else ''
	df = pandas.DataFrame.from_dict(result_dict).T

	df['rxn_flux'] = df['rxn_flux'].astype(float)
	df['met_flux'] = df['met_flux'].astype(float)

	df = df.loc[df['met_flux'].abs().sort_values(ascending=False).index]
	return df[df['ub'] != 0]

def get_reactions_of_met(me,met,s = 0, ignore_types = (),only_types = (), verbose = False,growth_key='mu'):
    met_stoich = 0
    if only_types:
        only_reaction_types = tuple([getattr(coralme.core.reaction,i) for i in only_types])
    elif ignore_types:
        ignore_reaction_types = tuple([getattr(coralme.core.reaction,i) for i in ignore_types])
    reactions = []

    if not hasattr(me.metabolites,met):
        return reactions
    for rxn in me.metabolites.get_by_id(met).reactions:
        if only_types and not isinstance(rxn, only_reaction_types):
            continue
        elif ignore_types and isinstance(rxn, ignore_reaction_types):
            continue
        try:
            met_obj = me.metabolites.get_by_id(met)
            pos = 1 if get_met_coeff(rxn.metabolites[met_obj],0.1,growth_key=growth_key) > 0 else -1
            rev = 1 if rxn.lower_bound < 0 else 0
            fwd = 1 if rxn.upper_bound > 0 else 0
        except:
            if verbose:
                print(rxn.id, ' could not parse')
            else:
                pass
        try:
            if not s:
                reactions.append(rxn)
                if verbose:
                    print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)

            elif s == pos*fwd or s == -pos*rev:
                reactions.append(rxn)
                if verbose:
                    print('(',rxn.id,rxn.lower_bound,rxn.upper_bound,')', '\t',rxn.reaction)

        except:
            if verbose:
                print(rxn.id, 'no reaction')
            else:
                pass
    return reactions

def find_issue(query,d,msg = ''):
    if isinstance(d,dict):
        if 'msg' in d:
            msg = d['msg']
            if 'triggered_by' in d:
                trigger = d['triggered_by']
                find_issue(query,trigger,msg=msg)
        else:
            for k,v in d.items():
                find_issue(query,v,msg=msg)
    elif isinstance(d,list):
        for i in d:
            find_issue(query,i,msg=msg)
    elif isinstance(d,str):
        if query == d:
            print(msg)
    else:
        raise TypeError("unsupported type  " + type(d))


# TODO: Add warnings
def fill_builder(b,fill_with='CPLX_dummy',key=None,d=None,fieldname=None,warnings=None):
    if isinstance(b,coralme.builder.main.MEBuilder):
        for i in dir(b.org):
            if i[0] == '_':
                continue
            attr = getattr(b.org,i)
            if not isinstance(attr,dict):
                continue
            fill_builder(attr,fill_with=fill_with,fieldname=i,warnings = warnings)
    elif isinstance(b,dict):
        for k,v in b.items():
            fill_builder(v,key=k,d=b,fill_with=fill_with,fieldname=fieldname,warnings=warnings)
    elif isinstance(b,list):
        include_keys = ['enzymes','proteins','enzyme','protein','machine']
        for ik in include_keys:
            if key in ik:
                if not b:
                    d[key] = ['CPLX_dummy']
    elif isinstance(b,str):
        include_keys = ['enzymes','proteins','enzyme','protein','machine']
        for ik in include_keys:
            if key in ik:
                if not b:
                    d[key] = 'CPLX_dummy'
    else:
        pass

def gap_find(me_model,de_type = None):
	#from draft_coralme.util.helper_functions import find_gaps

	logging.warning('  '*5 + 'Finding gaps in the ME-model...')
	me_gaps = coralme.builder.helper_functions.find_gaps(me_model, growth_key = me_model.mu)

	if de_type == 'me_only':
		logging.warning('  '*5 + 'Finding gaps from the M-model only...')
		m_gaps = coralme.builder.helper_functions.find_gaps(me_model.gem)
		idx = list(set(me_gaps.index) - set(m_gaps.index))
	else:
		idx = list(set(me_gaps.index))
	new_gaps = me_gaps.loc[idx]

	filt1 = new_gaps['p'] == 1
	filt2 = new_gaps['c'] == 1
	filt3 = new_gaps['u'] == 1

	deadends = list(new_gaps[filt1 | filt2 | filt3].index)
	deadends = sorted([ x for x in deadends if 'biomass' not in x if not x.endswith('_e') ])

	logging.warning('  '*5 + '{:d} metabolites were identified as deadends.'.format(len(deadends)))
	for met in deadends:
		name = me_model.metabolites.get_by_id(met).name
		logging.warning('  '*6 + '{:s}: {:s}'.format(met, 'Missing metabolite in the M-model.' if name == '' else name))
	return deadends

def gap_fill(me_model, deadends = [], growth_key_and_value = { sympy.Symbol('mu', positive = True) : 0.1 }, met_types = 'Metabolite'):
	if sys.platform in ['win32', 'darwin']:
		self.me_model.get_solution = self.me_model.optimize_windows
		self.me_model.get_feasibility = self.me_model.feas_gurobi
	else:
		self.me_model.get_solution = self.me_model.optimize
		self.me_model.get_feasibility = self.me_model.feasibility

	if len(deadends) != 0:
		logging.warning('  '*5 + 'Adding a sink reaction for each identified deadend metabolite...')
		coralme.builder.helper_functions.add_exchange_reactions(me_model, deadends, prefix='TS_')
	else:
		logging.warning('  '*5 + 'Empty set of deadends metabolites to test.')
		return None

	logging.warning('  '*5 + 'Optimizing gapfilled ME-model...')

	if me_model.get_feasibility(keys = growth_key_and_value):
		#logging.warning('  '*5 + 'The ME-model is feasible.')
		logging.warning('  '*5 + 'Gapfilled ME-model is feasible with growth rate {:g} 1/h.'.format(list(growth_key_and_value.values())[0]))
		return True
	else:
		#logging.warning('  '*5 + 'The ME-model is not feasible.')
		logging.warning('  '*5 + 'Provided set of sink reactions for deadend metabolites does not allow growth.')
		return False

def brute_force_check(me_model, metabolites_to_add, growth_key_and_value):
	if sys.platform == 'win32':
		me_model.get_solution = me_model.opt_gurobi
		me_model.get_feasibility = me_model.feas_gurobi
	else:
		me_model.get_solution = me_model.optimize
		me_model.get_feasibility = me_model.feasibility

	logging.warning('  '*5 + 'Adding sink reactions for {:d} metabolites...'.format(len(metabolites_to_add)))
# 	existing_sinks = [r.id for r in me_model.reactions.query('^TS_')]
	sk_rxns = coralme.builder.helper_functions.add_exchange_reactions(me_model, metabolites_to_add, prefix='TS_')

	if me_model.get_feasibility(keys = growth_key_and_value):
		pass
	else:
		logging.warning('  '*5 + 'Provided metabolites through sink reactions cannot recover growth. Proceeding to next set of metabolites.')
		return metabolites_to_add, [], False

	rxns = []
	rxns_to_drop = []
# 	rxns_to_append = []
# 	for idx, flux in me_model.solution.fluxes.items():
	for r in sk_rxns:
		idx = r.id
		flux = me_model.solution.fluxes[idx]
		if idx.startswith('TS_') and idx.split('TS_')[1] in metabolites_to_add:
# 			if r.id in existing_sinks:
# 				rxns_to_append.append(idx)
# 				continue
			if abs(flux) > 0:
				rxns.append(idx)
			else:
				#logging.warning('Closing {}'.format(idx))
				rxns_to_drop.append(idx)
				me_model.reactions.get_by_id(idx).bounds = (0, 0)

	logging.warning('  '*6 + 'Sink reactions shortlisted to {:d} metabolites.'.format(len(rxns)))

	# reaction_id:position in the model.reactions DictList object
# 	rxns = rxns + rxns_to_append# Try present SKs the last.
# 	logging.warning('  '*6 + 'Will try a total of {:d} metabolites including previous iterations:'.format(len(rxns)))
	ridx = []
	for r in rxns:
		ridx.append((r,me_model.reactions._dict[r]))
# 	ridx = { k:v for k,v in me_model.reactions._dict.items() if k in rxns }

	# populate with stoichiometry
	Sf, Se, lb, ub, b, c, cs, atoms, lambdas = me_model.construct_lp_problem()

	if lambdas is None:
		Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, Se, lb, ub, growth_key_and_value, atoms)
	else:
		Sf, Se, lb, ub = coralme.builder.helper_functions.evaluate_lp_problem(Sf, lambdas, lb, ub, growth_key_and_value, atoms)

	res = []
	msg = 'Processed: {:s}/{:d}, Gaps: {:d}. The ME-model is {:s}feasible if {:s} is closed.'
	for idx, (rxn, pos) in enumerate(ridx):
		lb[pos] = 0
		ub[pos] = 0
		if me_model.get_feasibility(keys = growth_key_and_value, **{'lp' : [Sf, dict(), lb, ub, b, c, cs, set(), lambdas]}):
			res.append(False)
			logging.warning('{:s} {:s}'.format('  '*6, msg.format(str(idx+1).rjust(len(str(len(ridx)))), len(ridx), len([ x for x in res if x ]), '', rxn)))
		else:
			lb[pos] = -1000
			ub[pos] = +1000
			res.append(True)
			logging.warning('{:s} {:s}'.format('  '*6, msg.format(str(idx+1).rjust(len(str(len(ridx)))), len(ridx), len([ x for x in res if x ]), 'not ', rxn)))

	bf_gaps = [ y for x,y in zip(res, rxns) if x ] # True
	no_gaps = [ y for x,y in zip(res, rxns) if not x ] + rxns_to_drop

	return bf_gaps, no_gaps, True

def get_mets_from_type(me_model,met_type):
	if met_type[1] == 'User guesses':
		return set(met_type[0])
	elif met_type == 'ME-Deadends':
		return set(coralme.builder.helper_functions.gap_find(me_model,de_type='me_only'))
	elif met_type == 'All-Deadends':
		return set(coralme.builder.helper_functions.gap_find(me_model))
	elif met_type == 'Cofactors':
		return set(coralme.builder.helper_functions.get_cofactors_in_me_model(me_model))
	else:
		mets = set()
		for met in me_model.metabolites:
			filter1 = type(met) == getattr(coralme.core.component, met_type)
			filter2 = met.id.startswith('trna')
			filter3 = met.id.endswith('trna_c')
			filter4 = met.id.endswith('_e')
			if filter1 and not filter2 and not filter3 and not filter4:
				mets.add(met.id)
		return mets

def _append_metabolites(mets,new_mets):
	return mets + [m for m in new_mets if m not in mets]

def brute_check(me_model, growth_key_and_value, met_type, skip = set(), history = dict()):
	mets = get_mets_from_type(me_model,met_type)
	if met_type == 'Metabolite':
		#remove from the metabolites to test that are fed into the model through transport reactions
		medium = set([ '{:s}_c'.format(x[3:-2]) for x in me_model.gem.medium.keys() ])
		mets = set(mets).difference(medium)
		# filter out manually
		mets = set(mets).difference(set(['ppi_c', 'ACP_c', 'h_c']))
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
	mets = set(mets).difference(skip)
	if met_type[1] == 'User guesses':
		history['User guesses'] = mets
	else:
		history[met_type] = mets

	mets_to_check = []
	for k,v in history.items():
		mets_to_check = _append_metabolites(mets_to_check,v)
	return history,coralme.builder.helper_functions.brute_force_check(me_model,
															  mets_to_check[::-1],
															  growth_key_and_value)

def evaluate_lp_problem(Sf, Se, lb, ub, keys, atoms):
	lb = [ x(*[ keys[x] for x in list(atoms) ]) if hasattr(x, '__call__') else float(x.xreplace(keys)) if hasattr(x, 'subs') else x for x in lb ]
	ub = [ x(*[ keys[x] for x in list(atoms) ]) if hasattr(x, '__call__') else float(x.xreplace(keys)) if hasattr(x, 'subs') else x for x in ub ]
	Se = { k:x(*[ keys[x] for x in list(atoms) ]) if hasattr(x, '__call__') else float(x.xreplace(keys)) if hasattr(x, 'subs') else x for k,x in Se.items() }

	Sf.update(Se)

	return Sf, Se, lb, ub

def substitute_value(m,coeff):
    if not hasattr(coeff,'subs'):
        return coeff
    return coeff.subs([(m._model.mu, 1e-6)])

def get_next_from_type(l,t):
    try:
        if isinstance(l,dict):
            return next(i for i,v in l.items() if isinstance(i,t) and substitute_value(i,v) > 0)
        return next(i for i in l if isinstance(i,t))
    except:
        return set()

def find_complexes(m, seen = set()):
    if not m:
        return set()
#     print(m.id)
    if m in seen:
        return set()
#     print(type(m))

    seen.add(m)

    # Metabolite objects
    if isinstance(m,coralme.core.component.TranslatedGene):
        cplxs = set()
        for r in m.reactions:
            if substitute_value(m,r.metabolites[m] > 0):
                continue
            cplxs = cplxs | find_complexes(r, seen=seen)
        return cplxs
    if isinstance(m,coralme.core.component.TranscribedGene):
        translated_protein = m.id.replace('RNA_','protein_')
        if translated_protein in m._model.metabolites:
            return find_complexes(m._model.metabolites.get_by_id(translated_protein), seen=seen)
        cplxs = set()
        for r in m.reactions:
            if substitute_value(m,r.metabolites[m] > 0):
                continue
            cplxs = cplxs | find_complexes(r, seen=seen)
        return cplxs
    if isinstance(m,coralme.core.component.ProcessedProtein):
        cplxs = set()
        for r in m.reactions:
            if substitute_value(m,r.metabolites[m] > 0):
                continue
            cplxs = cplxs | find_complexes(r, seen=seen)
        return cplxs

    if isinstance(m,coralme.core.component.Complex) or isinstance(m,coralme.core.component.GenericComponent) or isinstance(m,coralme.core.component.GenerictRNA):
        other_formations = [r for r in m.reactions if (isinstance(r,coralme.core.reaction.ComplexFormation) or isinstance(r,coralme.core.reaction.GenericFormationReaction)) and substitute_value(m,r.metabolites[m]) < 0]
#         print(other_formations)
        cplxs = set([m])
        if other_formations:
            for r in other_formations:
                cplxs = cplxs | find_complexes(r, seen=seen)
        return cplxs
#         return set([m])

    # Reaction objects
    if isinstance(m,coralme.core.reaction.PostTranslationReaction):
        return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.ProcessedProtein), seen=seen)
    if isinstance(m,coralme.core.reaction.ComplexFormation):
        return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen)
    if isinstance(m,coralme.core.reaction.GenericFormationReaction):
        return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen)
    if isinstance(m,coralme.core.reaction.tRNAChargingReaction):
        return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenerictRNA), seen=seen)
    if isinstance(m,coralme.core.reaction.MetabolicReaction):
        return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen) | \
                find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen)


    return set()

def get_functions(cplx):
	functions = set()
	for r in cplx.reactions:
		if isinstance(r,coralme.core.reaction.MetabolicReaction) and hasattr(r,'subsystem'):
			if r.subsystem:
				functions.add('Metabolic:' + r.subsystem)
				continue
			functions.add('Metabolic: No subsystem')
			continue
		if isinstance(r,coralme.core.reaction.TranslationReaction):
			functions.add('Translation')
		elif isinstance(r,coralme.core.reaction.TranscriptionReaction):
			functions.add('Transcription')
		elif isinstance(r,coralme.core.reaction.tRNAChargingReaction):
			functions.add('tRNA-Charging')
		elif isinstance(r,coralme.core.reaction.PostTranslationReaction):
			functions.add('Post-translation')
		elif isinstance(r,coralme.core.reaction.SummaryVariable):
			functions.add('Biomass')
	return functions

def dict_to_defaultdict(dct):
    return collections.defaultdict(lambda: [], dct)

def save_curation_notes(curation_notes,filepath):
	file = open(filepath,'w')
	file.write(json.dumps(curation_notes, indent=4))
	file.close()

def load_curation_notes(filepath):
	if not os.path.isfile(filepath):
		return collections.defaultdict(list)
	with open(filepath) as json_file:
		return json.load(json_file,object_hook=dict_to_defaultdict)

def publish_curation_notes(curation_notes,filepath):
	file = open(filepath,'w')
	for k,v in curation_notes.items():
		file.write('\n')
		for w in v:
			file.write('{} {}@{} {}\n'.format('#'*20,w['importance'],k,'#'*20))
			file.write('{} {}\n'.format('*'*10,w['msg']))
			if 'triggered_by' in w:
				file.write('The following items triggered the warning:\n')
				for i in w['triggered_by']:
					if isinstance(i,dict):
						file.write('\n')
						file.write(json.dumps(i))
						file.write('\n')
					else:
						file.write(i + '\n')
			file.write('\n{}Solution:\n{}\n\n'.format('*'*10,w['to_do']))
		file.write('\n\n')
	file.close()

def get_cofactors_in_me_model(me):
	cofactors = set()
	for i in me.process_data.query('^mod_'):
		for k,v in i.stoichiometry.items():
			if not me.metabolites.has_id(k):
				continue
			if v < 0:
				cofactors.add(k)
	return list(cofactors)

def is_producible(me,met,growth_key_and_value):
	if met not in me.metabolites:
		return False
	r = add_exchange_reactions(me, [met], prefix = 'DM_')[0]
	r.bounds = (1e-16,1e-16)
	if me.check_feasibility(keys = growth_key_and_value):
		r.remove_from_model()
		return True
	else:
		r.remove_from_model()
		return False


def get_transport_reactions(model,met_id,comps=['e','c']):
    from_met = re.sub('_[a-z]$','_'+comps[0],met_id)
    to_met = re.sub('_[a-z]$','_'+comps[1],met_id)

    if isinstance(model,coralme.core.model.MEModel):
        reaction_type = ['MetabolicReaction']
    else:
        reaction_type = 0
    prod_rxns = [rxn.id for rxn in get_reactions_of_met(model,to_met,s=1,verbose=0,only_types=reaction_type)]
    cons_rxns = [rxn.id for rxn in get_reactions_of_met(model,from_met,s=-1,verbose=0,only_types=reaction_type)]
    transport_rxn_ids = list(set(prod_rxns)&set(cons_rxns))

    return [model.reactions.get_by_id(rxn_id) for rxn_id in transport_rxn_ids]

def get_all_transport_of_model(model):
    transport_reactions = []
    for r in tqdm.tqdm(model.reactions):
        comps = r.get_compartments()
        if len(comps) > 1:
            transport_reactions.append(r.id)
    return list(set(transport_reactions))


def format_kcats_from_DLKcat(df):
# df = pandas.read_csv("./bacillus/building_data/bacillus_rxn_kcats.tsv",sep='\t',index_col=0).set_index("reaction")

    df2 = pandas.DataFrame(columns=["direction","complex","mods","keff"])
    for r,keff in df['Kcat value (1/s)'].items():
        d = {}
        rid,cplx = re.split("_FWD_|_REV_",r)
        base_cplx = cplx.split("_mod_")[0]
        mods = cplx.split(base_cplx)[1]
        mods = " AND ".join(mods.split("_mod_")[1:])
        direc = "FWD" if "FWD" in r else "REV"
        d[rid] = {
            "direction":direc,
            "complex":base_cplx,
            "mods":mods,
            "keff":keff
        }
        df2 = pandas.concat([df2,pandas.DataFrame.from_dict(d).T])
    df2.index.name = 'reaction'
    return df2
#     df2.to_csv("./bacillus/building_data/keffs.txt",sep='\t')
