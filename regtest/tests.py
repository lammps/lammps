"""

LAMMPS regression testing package.
Author: Peter Shannon, pietro.shannon@gmail.com

Module contains tests.

"""
import re
import os
import sys
import math

def single_input(file, hook):
	"""
	Returns lines from file between hooks.
	"""
	
	input_hook = []

	hook_begin = False
	hook_end = False

	try:

		single_file = open(file)

		for line in single_file:

				if ("#" + hook) in line and hook_begin is True:
					hook_end = True

				if hook_begin is True and hook_end is False:
					input_hook.append(line)

				if ("#" + hook) in line and hook_begin is False:
					hook_begin = True

		single_file.close()
		
		return input_hook

	except IOError:
		
		print "The file " + file + " does not exist."

#Extracts numbers from LAMMPS run output
def run_numbers(lists):
	"""
	Returns a list which contains the numbers generated from 
	the "run" command.  This configuration currently does not
	support "thermo_style one" and "thermo_style multi" commands. 

	"""
	
	numbers = []

	number_begin = False
	number_end = False

	for i in lists:

		if "Loop" in i and number_begin is True:
			number_end = True

		if number_begin is True and number_end is False:
			numbers.append(i.split())
	
		if "Step" in i and number_begin is False:
			number_begin = True

	return numbers

def test_run_numbers(default, test, label, rcoef, filename):
	"""
	Prints results of a compairson between two files with numbers
	generated from the "run" command.  
	"""

	status = None
	errors = []
	
	for i in range(len(default)):
		for j in range(len(default[i])):
			coef = math.fabs(float(default[i][j]) * float(rcoef[j]))
			if float(test[i][j]) < float(default[i][j]) - coef or float(test[i][j]) > float(default[i][j]) + coef: 
				errors.append(str(label[j]) + " mismatch: " + str(test[i][j]) + " instead of " + str(default[i][j]) + \
				" + or - " + str(coef) + " at step " + str(default[i][0]))
				status = "Fail"
	if status is not "Fail":
		status = "Pass"

	print filename +": " + status  
	for e in errors: 
		print "  -" + e


