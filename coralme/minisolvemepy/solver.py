#!/usr/bin/python3

import numpy
import scipy

import cobra
from cobra.core import Solution
import coralme
from coralme.minisolvemepy import qwarmLP, warmLP

# from solvemepy.me1
from sympy import lambdify, Basic, Symbol

# Modified from solvemepy.me2
class ME_NLP:
    """
    Contains the data matrices needed for solving ME as an NLP using qMINOS
    """
    def __init__(self, me, growth_rxn = 'biomass_dilution'):
        # The ME model object
        self.me = me

        # Reformulation of ME to NLP
        self.A = None
        self.B = None
        self.S = None
        self.b = None
        self.c = None
        self.xl= None
        self.xu= None
        # Inputs to qminos
        self.J = None
        self.nnCon = None
        self.nnJac = None
        self.neJac = None
        self.ne    = None
        self.ha    = None
        self.ka    = None
        self.ad    = None
        self.bld   = None
        self.bud   = None
        self.mu0   = None
        self.probname = "me_nlp"
        self.M        = None
        self.N        = None
        self.nb       = None
        # Solution and exit flag
        self.x      = None
        self.inform = numpy.array(0)
        # Hold LP results and options
        self.lp_inform  = None
        self.lp_hs      = None
        self.lp_x       = None
        # Initialize solver options
        self.init_solver_opts()

        # from solvemepy.me1
        self.substitution_dict = { self.me.mu: 1. }
        # Initially have None for compiled_expressions
        self.compiled_expressions = None

    def init_solver_opts(self):
        #----------------------------------------------------
        # Solver options
        self.opt_strwidth  = {}
        self.opt_realwidth = {}
        self.opt_intwidth  = {}
        self.opt_strlist   = {}
        self.opt_intdict   = {}
        self.opt_realdict  = {}
        self.opt_stropts   = {}
        self.opt_intopts   = {}
        self.opt_intvals   = {}
        self.opt_realopts  = {}
        self.opt_realvals  = {}
        #----------------------------------------------------
        # NLP solver options
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['nlp'] = 72
        self.opt_realwidth['nlp'] = 55
        self.opt_intwidth['nlp'] = 55
        self.opt_strlist['nlp'] = [
                'Maximize',     # Default obj sense is to maximize
                'Completion full',
                'Print level (jflxb) 00001',
                'Solution No'
                ]
        self.opt_intdict['nlp'] = {
                'Major iterations': 1000,
                'Superbasics limit': 40,
                'Verify level': 0,
                'Scale option': 2,
                'Partial price': 1,
                'Iterations': 10000,
                'Print frequency': 100000,
                'Summary level': 0,
                'Summary frequency': 100,
                'Solution file': 9,
                'New basis file': 11,
                'Save frequency': 500000
                }
        self.opt_realdict['nlp'] = {
                'Penalty parameter':100.0,
                'LU factor tol': 1.1,
                'LU update tol': 1.1,
                'LU singularity tol': 1e-30,
                'Feasibility tol': 1e-15,
                'Optimality tol': 1e-15,
                'Unbounded step size': 1e+30
                }

        #----------------------------------------------------
        # LP options
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['lp'] = 72
        self.opt_realwidth['lp'] = 55
        self.opt_intwidth['lp'] = 55
        self.opt_strlist['lp'] = [
                'Maximize',     # Default obj sense is to maximize
                'Solution No'
                ]
        self.opt_intdict['lp'] = {
                'New basis file': 11,
                'Save frequency': 500000,
                'Print level': 0,
                'Print frequency': 100000,
                'Scale option': 2,
                'Iteration limit': 2000000,
                'Expand frequency': 100000
                }
        self.opt_realdict['lp'] = {
                'Penalty parameter':100.0,
                'LU factor tol': 10.0,
                'LU update tol': 10.0,
                'LU singularity tol': 1e-30,
                'Feasibility tol': 1e-20,
                'Optimality tol': 1e-20,
                'Unbounded step size': 1e+30
                }

        #----------------------------------------------------
        # LP options: double-precision (can't set as strict tols)
        #----------------------------------------------------
        # Width of characters allowed in each options for qMINOS
        self.opt_strwidth['lp_d'] = 72
        self.opt_realwidth['lp_d'] = 55
        self.opt_intwidth['lp_d'] = 55
        self.opt_strlist['lp_d'] = [
                'Maximize',     # Default obj sense is to maximize
                'Solution No'
                ]
        self.opt_intdict['lp_d'] = {
                'New basis file': 11,
                'Save frequency': 500000,
                'Print level': 0,
                'Print frequency': 100000,
                'Scale option': 2,
                'Iteration limit': 2000000,
                'Expand frequency': 100000
                }
        self.opt_realdict['lp_d'] = {
                'Penalty parameter':100.0,
                'LU factor tol': 1.9,
                'LU update tol': 1.9,
                'LU singularity tol': 1e-12,
                'Feasibility tol': 1e-7,
                'Optimality tol': 1e-7,
                'Unbounded step size': 1e+18
                }

    def get_solver_opts(self, prob = 'lp'):
        """
        Return options that will be passed as arguments to minoss
        """

        lst = [ c for c in [ s.ljust(self.opt_strwidth[prob]) for s in self.opt_strlist[prob] ] ]
        stropts = numpy.array(numpy.array(lst, dtype = 'c').T)

        intkeys = self.opt_intdict[prob].keys()
        realkeys = self.opt_realdict[prob].keys()

        lst = [ c for c in [ s.ljust(self.opt_intwidth[prob]) for s in intkeys ] ]
        intopts = numpy.array(numpy.array(lst, dtype = 'c').T)

        lst = [ c for c in [ s.ljust(self.opt_realwidth[prob]) for s in realkeys ] ]
        realopts = numpy.array(numpy.array(lst, dtype = 'c').T)

        intvals = numpy.array([ self.opt_intdict[prob][k] for k in intkeys ], dtype = 'i4')
        realvals = numpy.array([ self.opt_realdict[prob][k] for k in realkeys ], dtype = 'd')

        self.opt_stropts[prob] = stropts
        self.opt_intopts[prob] = intopts
        self.opt_realopts[prob] = realopts
        self.opt_intvals[prob] = intvals
        self.opt_realvals[prob]= realvals

        nStrOpts = len(self.opt_strlist[prob])
        nIntOpts = len(self.opt_intdict[prob].keys())
        nRealOpts = len(self.opt_realdict[prob].keys())

        return stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts

    def makeME_LP(self, S, b, c, xl, xu, csense):
        """
        Create simple LP for qMINOS and MINOS
        Inputs:
        nlp_compat  Make matrices compatible with NLP so that basis can
                    be used to warm start NLP by setting
        12 Aug 2015: first version
        """

        # c is added as a free (unbounded slacks) row,
        # so that MINOS treats problem as an LP - Ding Ma
        J = scipy.sparse.vstack((S, c)).tocsc()
        J.sort_indices()

        if hasattr(b, 'tolist'):
            b = b.tolist()

        b2 = b + [0.0]
        m, n = J.shape
        ne = J.nnz
        # Finally, make the P, I, J, V, as well
        # Row indices: recall fortran is 1-based indexing
        I = [ i+1 for i in J.indices ]
        V = J.data
        # Pointers to start of each column
        # Just change to 1-based indexing for Fortran
        P = [ pi+1 for pi in J.indptr ]

        # Make primal and slack bounds
        bigbnd = 1e+40
        # For csense==E rows (equality)
        sl = numpy.matrix([ bi for bi in b2 ]).transpose()
        su = numpy.matrix([ bi for bi in b2 ]).transpose()

        for row, csen in enumerate(csense):
            if csen == 'L':
                sl[row] = -bigbnd
            elif csen == 'G':
                su[row] = +bigbnd

        # Objective row has free bounds
        sl[m-1] = -bigbnd
        su[m-1] = +bigbnd

        bl = scipy.vstack([xl, sl])
        bu = scipy.vstack([xu, su])

        return J, ne, P, I, V, bl, bu

    # from solvemepy.me2
    def construct_S(self, growth_rate):
        """
        From cobrame--in case me does not have construct_S
        """
        me = self.me
        growth_key = self.me.mu

        # intialize to 0
        S = scipy.sparse.dok_matrix((len(me.metabolites), len(me.reactions)))
        # populate with stoichiometry
        for i, r in enumerate(me.reactions):
            for met, value in r._metabolites.items():
                met_index = me.metabolites.index(met)
                if hasattr(value, "subs"):
                    S[met_index, i] = float(value.subs(growth_key, growth_rate))
                else:
                    S[met_index, i] = float(value)
        return S

    # from solvemepy.me1
    def compile_expr(self, expr):
        """Compile expressions with all parameter keys in ME 1.0"""
        # Note that ufuncify too slow to compile
        #f = lambdify(self.subs_keys_ordered, expr) if isinstance(expr, Basic) else expr
        # 19 Jan 2017:  already checked isinstance(expr, Basic) before calling this method
        f = lambdify(self.subs_keys_ordered, expr)
        return f

    def compile_expressions(self):
        """
        Compile expressions for ME 1.0.
        Use format consistent with cobrame:
        (met_index, rxn_index): stoichiometry, (None, rxn_index): (lower_bound, upper_bound)
        (met_index, None): (met_bound, met_constraint_sense)
        """

        expressions = {}
        me = self.me

        # Reaction bounds
        for idx, rxn in enumerate(me.reactions):
            for met, stoic in rxn._metabolites.items():
                if isinstance(stoic, Basic):
                    expressions[(me.metabolites.index(met), idx)] = self.compile_expr(stoic)

            # If lower or upper bound symbolic:
            if isinstance(rxn.lower_bound, Basic) or isinstance(rxn.upper_bound, Basic):
                expressions[(None, idx)] = (self.compile_expr(rxn.lower_bound), self.compile_expr(rxn.upper_bound))

        # Metabolite bound
        for idx, metabolite in enumerate(me.metabolites):
            if isinstance(metabolite._bound, Basic):
                expressions[(idx, None)] = (self.compile_expr(metabolite._bound), metabolite._constraint_sense)

        return expressions

    def make_lp(self, mu_fix, compiled_expressions = None):
        """
        Construct LP problem for qMINOS or MINOS
        """

        me = self.me
        if self.compiled_expressions is None:
            # 16 Sep 2016: [LY] should update ordered keys, too, in case the reason
            # compiled_expressions is None is because new symbols were added
            self.subs_keys_ordered = self.substitution_dict.keys()
            self.compiled_expressions = self.compile_expressions()

        self.substitution_dict[self.me.mu] = mu_fix
        # Substitution keys need to be the same as when the lambdas were originally compiled

        # Get the subtitution values in the right order for lambdify
        sub_vals = [self.substitution_dict[k] for k in self.subs_keys_ordered]

        # Initialize S, lb, ub, b
        S = scipy.sparse.dok_matrix((len(me.metabolites), len(me.reactions)))
        xl = numpy.matrix([r.lower_bound for r in me.reactions]).transpose()
        xu = numpy.matrix([r.upper_bound for r in me.reactions]).transpose()
        b = [0. for m in me.metabolites]

        # Fill in all matrix & constraint rhs entries (incl. not mu-dependent)
        for mind,met in enumerate(me.metabolites):
            # Fill in constraint bounds: MOSTLY just float
            if hasattr(met._bound, 'subs'):
                expr = self.compiled_expressions[(mind,None)]
                b[mind] = float(expr[0](*sub_vals))
            else:
                b[mind] = met._bound

            # Fill in stoichiometries: MOSTLY symbolic, or float? Hard to say.
            for rxn in met.reactions:
                rind = me.reactions.index(rxn)
                if (mind, rind) in self.compiled_expressions:
                    expr = self.compiled_expressions[(mind,rind)]
                    #****************************************
                    # DEBUG
                    try:
                        s = float(expr(*sub_vals))
                    except TypeError as e:
                        # Just indicate which rxn,met,stoich had issues
                        print(repr(e))
                        print('rxn=%s \t met=%s' % (rxn.id, met.id))
                        print('stoich=', rxn.metabolites[met])
                        raise Exception('Failed to convert symbolic stoichiometry to float')

                    if not numpy.isinf(s):
                        S[mind, rind] = s
                else:
                    S[mind,rind] = rxn.metabolites[met]

        # Fill in var bounds: MOSTLY just float
        for rind,rxn in enumerate(me.reactions):
            if hasattr(rxn.lower_bound, 'subs'):
                # Then, there must be a compiled expression
                expr = self.compiled_expressions[(None,rind)]
                xl[rind] = float(expr[0](*sub_vals))
            else:
                xl[rind] = rxn.lower_bound

            if hasattr(rxn.upper_bound, 'subs'):
                expr = self.compiled_expressions[(None,rind)]
                xu[rind] = float(expr[1](*sub_vals))
            else:
                xu[rind] = rxn.upper_bound

        c = [r.objective_coefficient for r in me.reactions]
        csense = ['E' for m in me.metabolites]

        J, ne, P, I, V, bl, bu = self.makeME_LP(S, b, c, xl, xu, csense)

        # Solve a single LP
        m,n = J.shape
        ha = I
        ka = P
        ad = V
        bld = [bi for bi in bl.flat]
        bud = [bi for bi in bu.flat]
        nb = m + n
        hs = numpy.zeros(nb, numpy.dtype('i4'))
        return m, n, ha, ka, ad, bld, bud, hs

    # from solvemepy.me2
    #def make_lp(self, mu_fix):
        #"""
        #Construct LP problem for qMINOS or MINOS.
        #"""

        #me = self.me
        #S = self.construct_S(mu_fix).tocsc()
        ##S = me.construct_s_matrix(mu_fix).tocsc()
        #xl = numpy.matrix([ r.lower_bound for r in me.reactions ]).transpose()
        #xu = numpy.matrix([ r.upper_bound for r in me.reactions ]).transpose()

        ## Also substitute mu in bounds
        #for idx, rxn in enumerate(me.reactions):
            #lb = rxn.lower_bound
            #ub = rxn.upper_bound
            #if hasattr(lb, 'subs'):
                ##xl[idx] = float(lb.subs(self.mu, mu_fix))
                #xl[idx] = float(lb.subs(me.mu, mu_fix))
            #if hasattr(ub, 'subs'):
                ##xu[idx] = float(ub.subs(self.mu, mu_fix))
                #xu[idx] = float(ub.subs(me.mu, mu_fix))

        ## S * v = b
        #b = [ m._bound for m in me.metabolites ]
        #c = [ r.objective_coefficient for r in me.reactions ]
        ## constraint sense eventually be in the metabolite...
        #csense = [ 'E' for m in me.metabolites ]
        ##csense = [m._constraint_sense for m in me.metabolites]
        #J, ne, P, I, V, bl, bu = self.makeME_LP(S, b, c, xl, xu, csense)

        ## Solve a single LP
        #m, n = J.shape
        #ha = I
        #ka = P
        #ad = V
        #bld = [ bi for bi in bl.flat ]
        #bud = [ bi for bi in bu.flat ]
        #nb = m + n
        #hs = numpy.zeros(nb, numpy.dtype('i4'))
        #return m, n, ha, ka, ad, bld, bud, hs

    def solvelp(self, muf, basis = None, precision = 'quad'):
        """
        x, status, hs = solvelp(self, muf, basis = None, precision = 'quad')

        Solve LP at mu using qMINOS or MINOS.
        Pass the basis (hs) back and forth with Fortran for warm-start.

        Inputs:
        muf: fixed growth rate
        basis: basis vector

        Outputs:
        x: primal solution
        status: solver status
        hs: basis
        """
        me = self.me

        m, n, ha, ka, ad, bld, bud, hs0 = self.make_lp(muf)

        hs = basis
        if hs is None:
            warm = False
            hs = hs0
        else:
            warm = True

        inform = numpy.array(0)
        probname = 'me_lp'
        precision = precision.lower()

        if precision == 'quad':
            optimizer = qwarmLP.qwarmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp')

        elif precision == 'double':
            optimizer = warmLP.warmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp_d')

        elif precision == 'dq' or precision == 'dqq':
            # D
            self.opt_intdict['lp_d']['Scale option'] = 2
            optimizer = warmLP.warmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp_d')

            # Q1: pass optimal basis hs and scale = 2
            warm = True
            self.opt_intdict['lp']['Scale option'] = 2
            optimizer = qwarmLP.qwarmlp
            stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp')

            # Last Q2 if requested: pass optimal basis hs and scale = 0
            if precision == 'dqq':
                self.opt_intdict['lp']['Scale option'] = 0
                optimizer = qwarmLP.qwarmlp
                stropts, intopts, realopts, intvals, realvals, nStrOpts, nIntOpts, nRealOpts = self.get_solver_opts('lp')

                # Kindly reset scale option to default
                self.opt_intdict['lp']['Scale option'] = 2

        else:
            raise ValueError('The \'precision\' must be \'quad\', \'double\', \'dq\', or \'dqq\'. Provided: {:s}'.format(str(precision)))

        x, pi, rc = optimizer(
            inform, probname, m, ha, ka, ad, bld, bud, hs, warm,
            stropts, intopts, realopts, intvals, realvals,
            nstropts = nStrOpts, nintopts = nIntOpts, nrealopts = nRealOpts
            )

        self.inform = inform
        self.hs = hs
        self.lp_hs = hs
        self.x = x
        # Save dual and reduced cost information
        self.pi = pi
        # Reduced cost: g - (A I)'*pi, where g is the gradient, l <= A*x <= u are constraints
        # including the objective function in the last row
        self.rc = rc

        # Write the solution to the ME model's solution for a consistent solve interface
        #f = x[0]
        # Aug 27, 2015: obj coeffs are not always mu (or x[0])
        f = sum([ rxn.objective_coefficient * x[idx] for idx, rxn in enumerate(self.me.reactions) ])
        x_primal = x[0:len(self.me.reactions)]   # The remainder are the slacks
        x_dict = { rxn.id:x[idx] for idx, rxn in enumerate(self.me.reactions) }
        y = pi
        # J = [S; c]
        y_dict = { met.id:y[idx] for idx, met in enumerate(me.metabolites) }
        #y_dict['linear_objective'] = y[len(y)-1]

        status = self.inform
        if int(status) == 0:
            status = 'optimal'
        #self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)
        self.me.solution = cobra.core.Solution(
            objective_value = x_dict['biomass_dilution'],
            status = status,
            fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
            reduced_costs = y_dict,
            shadow_prices = None,
            )

        return x, status, hs

    def bisectmu(
        self, precision = 1e-6, mumin = 0.0, mumax = 2.0, maxIter = 100,
        basis = None, verbosity = False, solver_precision = 'quad'
        ):
        """
        muopt, hs, xopt, cache = bisectmu(
            self, precision = 1e-6, mumin = 0.0, mumax = 2.0, maxIter = 100,
            basis = None, verbosity = False, solver_precision = 'quad'
            )

        Bisection to maximize mu using qMINOS.
        Sequence of feasibility problems.
        """
        import copy as cp

        me = self.me
        hs = basis

        # basis (hs) intent(inout). Will generate first basis from Cold-start if hs=None
        # Save solutions, recovered at golden ratios
        cache = {}
        # Store the final Solution object
        solution = None

        def checkmu(muf, hs):
            if muf not in cache:
                x_new, stat_new, hs_new = self.solvelp(muf, basis = hs, precision = solver_precision)
                if me.solution.status == 'optimal':
                    hs = hs_new
                stat = me.solution.status
                sol = cp.deepcopy(me.solution)
                cache[muf] = stat

            return cache[muf], hs, sol, x_new

        muopt = mumin
        xopt = None

        if verbosity:
            print('Iteration\t       Growth Rate\t Solution to check\tSolver Status')
            print('---------\t------------------\t------------------\t-------------')

        for idx in range(1, maxIter + 1):
            # Just a sequence of feasibility checks
            mu1 = (mumin + mumax) / 2.
            # Retrieve evaluation from cache if it exists: golden section advantage
            stat1, hs, sol1, x1 = checkmu(mu1, hs)
            if stat1 == 'optimal':
                mumin = mu1
                muopt = mu1
                solution = sol1
                xopt = x1
            else:
                mumax = mu1

            if verbosity:
                print('{:s}\t{:.16f}\t{:.16f}\t{:s}'.format(
                    str(idx).rjust(9), muopt, mu1, 'Not feasible' if stat1 == 1 else stat1.capitalize()))

            if abs(mumax - mumin) <= precision:
                break

        # Save final solution
        me.solution = solution
        # Save feasible basis
        self.feas_basis = hs

        return muopt, hs, xopt, cache
