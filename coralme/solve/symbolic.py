import coralme
import sympy

def _compile(expr, variable = coralme.util.mu):
	"""Compiles a sympy expression"""
	return sympy.lambdify(variable, expr) if isinstance(expr, sympy.Basic) else expr

def _eval(expr, mu0):
	"""
	Evaluate the expression as a function of mu
	will return the expr if it is not a function
	"""
	return expr(mu0) if callable(expr) else expr

def compile_expressions(me_model, variable = coralme.util.mu):
	"""
	Compiles symbolic expressions of mu to functions.
	The compiled expressions dict has the following key value pairs:
	(met_index, rxn_index): stoichiometry,
	(None, rxn_index): (lower_bound, upper_bound)
	(met_index, None): (met_bound, met_constraint_sense)
	"""
	expressions = {}

	# Reactions
	for i, r in enumerate(me_model.reactions):
		# stoichiometry
		#for met, stoic in iteritems(r._metabolites):
		for met, stoic in r._metabolites.items():
			if isinstance(stoic, sympy.Basic):
				expressions[(me_model.metabolites.index(met), i)] = sympy.lambdify(variable, stoic)
		# If either the lower or upper reaction bounds are symbolic
		if isinstance(r.lower_bound, sympy.Basic) or isinstance(r.upper_bound, sympy.Basic):
			expressions[(None, i)] = (_compile(r.lower_bound, variable), _compile(r.upper_bound, variable))

	# Metabolite bound
	for i, metabolite in enumerate(me_model.metabolites):
		if isinstance(metabolite._bound, sympy.Basic):
			expressions[(i, None)] = (_compile(metabolite._bound, variable), metabolite._constraint_sense)

	return expressions

def substitute_mu(lp, mu_value, compiled_expressions, solver_module = None):
	"""
	Substitute mu into a constructed LP
	mu: float
	"""
	# This only works for object-oriented solver interfaces. For other
	# solvers, need to pass in solver_module
	if solver_module is None:
		solver_module = lp.__class__
	#for index, expr in iteritems(compiled_expressions):
	for index, expr in compiled_expressions.items():
		# reaction bounds
		if index[0] is None:
			solver_module.change_variable_bounds(lp, index[1], _eval(expr[0], mu_value), _eval(expr[1], mu_value))

		# metabolite _bound
		elif index[1] is None:
			solver_module.change_constraint(lp, index[0], expr[1], _eval(expr[0], mu_value))

		# stoichiometry
		else:
			solver_module.change_coefficient(lp, index[0], index[1], _eval(expr, mu_value))
