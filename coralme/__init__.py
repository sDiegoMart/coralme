#!/usr/bin/python3
__version__ = "1.0"

import coralme.builder
import coralme.core
import coralme.io

import sys
if sys.platform in ['win32', 'darwin']:
	pass
else:
	import coralme.solver.solver

import coralme.util
