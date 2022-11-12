#!/usr/bin/python3
__version__ = "1.0"

# If cobrapy and optlang cannot handle symbolic parameters, assume using cobrapy versions <= 0.5.11
try:
	from optlang.interface import SymbolicParameter
except ImportError:
	import sympy
	mu = sympy.Symbol('mu', positive = True)
else:
	mu = SymbolicParameter('mu', value = 1.)
	mu._assumptions._tell('positive', True)
	mu._assumptions._tell('nonzero', True)
	mu._assumptions._tell('negative', False)
	mu._assumptions._tell('zero', False)
	mu._assumptions._tell('real', True)
	# mu._assumptions._tell('nonnegative', True)
	mu._assumptions['uuid'] = None

import coralme.util.building
import coralme.util.dogma
import coralme.util.massbalance
import coralme.util.me_model_interface
