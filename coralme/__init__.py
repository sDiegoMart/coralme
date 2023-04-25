#!/usr/bin/python3
__version__ = "1.0"

import coralme.builder
import coralme.core
import coralme.io

import sys
if sys.platform == 'win32':
	pass
else:
	import coralme.solver.solver

import coralme.util
