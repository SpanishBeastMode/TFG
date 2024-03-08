from pysmt.shortcuts import *
from pysmt.rewritings import nnf

index = 0

def tseitin(formula):
	global index
	if formula.is_symbol(BOOL) or formula.is_not():
		#print("SYMBOL", formula)
		return formula, []
	elif formula.is_and():
		t = Symbol("_T{}".format(index), BOOL)
		index += 1
		#print("CHECK AND FORMULA", formula)
		t1, l = tseitin(formula.arg(0))
		t2, r = tseitin(formula.arg(1))
		return t, l + r + [Or(Not(t), t1), Or(Not(t), t2), Or(t, Not(t1), Not(t2))]
		# return t, And(Or(Not(t), t1), Or(Not(t), t2), Or(t, Not(t1), Not(t2)),l,r)
	elif formula.is_or():
		t = Symbol("_T{}".format(index), BOOL)
		index += 1
		#print("CHECK OR FORMULA", formula)
		t1, l = tseitin(formula.arg(0))
		t2, r = tseitin(formula.arg(1))
		return t, l + r + [Or(Not(t), t1, t2), Or(t, Not(t1)), Or(t, Not(t2))]
		#return t, And(Or(Not(t), t1, t2), Or(t, Not(t1)), Or(t, Not(t2)),l,r)

def _plra_rec(formula, pos_polarity):
	"""This method extract all sub formulas in the formula and returns them as a dictionary.

	Args:
		formula (FNode): The formula to parse.
		pos_polarity (bool): The polarity of the formula.

	Returns:
		dict: the list of FNode in the formula with the corresponding truth value.
		bool: boolean that indicates if there are no more truth assignment to extract.

	"""
	over = True
	ass = {}
	if formula.is_or():
		return {}, True
	for arg in formula.args():
		#print(arg)
		if arg.is_symbol():
			over = False
			#print("NOT OVER")
			ass[arg] = Bool(True) 
		elif arg.is_not():
			over = False
			#print("NOT OVER")
			ass[arg.arg(0)] = Bool(False)
	return ass, over


def unit_propagation(formula, inputs):
	# A = Symbol("A", BOOL)
	# B = Symbol("B", BOOL)
	# C = Symbol("C", BOOL)
	# D = Symbol("D", BOOL)
	# E = Symbol("E", BOOL)
	# formula = And(And(C, Iff(Not(A),C)), Or(D, E))
	# inputs={E: False}
	# input_assignments = {variable: Bool(value) for variable, value in inputs.items()}

	formula = nnf(formula)
	# print(formula)

	t0, formula = tseitin(formula)

	# print(formula)

	formula = simplify(And(formula))

	# print(formula)
	# print("**************")

	f_before = formula
	f_next = formula
	atom_assignments = {t0: Bool(True), **inputs}
	lra_assignments = atom_assignments

	# iteratively simplify F[A<-mu^A], getting (possibily part.) mu^LRA
	while True:
		f_before = f_next
		f_next = simplify(substitute(f_before, lra_assignments))
		#print(f_next)
		lra_assignments, over = _plra_rec(f_next, True)
		#print(lra_assignments, over)
		# subs = {k: Bool(v) for k, v in lra_assignments.items()}
		atom_assignments.update(lra_assignments)
		#print("atom_assignments", atom_assignments)
		if over or lra_assignments == {}:
			break
	atom_assignments = {x:atom_assignments[x] for x in atom_assignments if x.symbol_name()[0] != "_"}
	# print(atom_assignments)
	return atom_assignments