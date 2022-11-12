"""
Optlang Interface for ME-model. Overrides some of the functions in the
high-level interface. If version of cobrapy < 0.6.0, this import does nothing.
"""

import six

try: from optlang.interface import SymbolicParameter
except ImportError:
	pass
else:
	# cobra/util/context.py:111
	# cobra/core/model.py:143
	# cobra/util/solver.py:348
	# <module 'coralme.util.me_model_interface' from
	# 'coralme/util/me_model_interface.py'> is not a valid solver interface. Pick one
	# from glpk_exact, glpk, gurobi, scipy.
	import optlang

	@six.add_metaclass(optlang.util.inheritdocstring)
	class Variable(optlang.interface.Variable):
		def __init__(self, name, lb = None, ub = None, type = 'continuous', *args, **kwargs):
			if type != 'continuous':
				raise ValueError('ME-models require continuous variables.')
			super(Variable, self).__init__(name, lb, ub, type, *args, **kwargs)

	@six.add_metaclass(optlang.util.inheritdocstring)
	class Constraint(optlang.interface.Constraint):
		def __init__(self, expression, sloppy = False, *args, **kwargs):
			super(Constraint, self).__init__(expression, sloppy = sloppy, *args, **kwargs)
		def set_linear_coefficients(self, coefficients):
			return

	@six.add_metaclass(optlang.util.inheritdocstring)
	class Objective(optlang.interface.Objective):
		def __init__(self, expression, sloppy = False, **kwargs):
			super(Objective, self).__init__(expression, sloppy = sloppy, **kwargs)
		def set_linear_coefficients(self, coefficients):
			return

	@six.add_metaclass(optlang.util.inheritdocstring)
	class OptimizationExpression(optlang.interface.OptimizationExpression):
		def __init__(self, expression, sloppy = False, **kwargs):
			super(OptimizationExpression, self).__init__(expression, sloppy = sloppy, **kwargs)
		def set_linear_coefficients(self, coefficients):
			return

	@six.add_metaclass(optlang.util.inheritdocstring)
	class Configuration(optlang.interface.MathematicalProgrammingConfiguration):
		def __init__(self, *args, **kwargs):
			super(Configuration, self).__init__(*args, **kwargs)

	@six.add_metaclass(optlang.util.inheritdocstring)
	class Model(optlang.interface.Model):
		def __init__(self, problem = None, *args, **kwargs):
			super(Model, self).__init__(*args, **kwargs)
