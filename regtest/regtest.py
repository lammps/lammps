"""

LAMMPS regression testing package.
Author: Peter Shannon, pietro.shannon@gmail.com

Run regression tests for lammps.

"""
import re
import tests
import os
import sys

class Config():

	default = None
	test = None
	test_folder = None
	hooks = []
	run_data = []
	run_coef = []

	def __init__(self):
		
		try:

			config = open("config")
			for line in config:
				
				if re.match("default =", line):
					self.default = line.split()[2]

				if re.match("test =", line):
					self.test = line.split()[2]

				if re.match("test_folder =", line):
					self.test_folder = line.split()[2]

				if re.match("hooks =", line):
					words = line.split()
					for i in range(2, len(words)):
						self.hooks.append(words[i])

				if re.match("run_data =", line):
					words = line.split()
					for i in range(2, len(words)):
						self.run_data.append(words[i])

				if re.match("run_coef =", line):
					words = line.split()
					for i in range(2, len(words)):
						self.run_coef.append(words[i])

			config.close()

		except IOError:
			print "No config file found in root directory"
			sys.exit()

"""
#############################################
					Main
#############################################
"""

config = Config()

default = tests.run_numbers(tests.single_input(config.default, config.hooks[0]))

if config.test_folder is not None:
	
	for files in os.listdir(str(config.test_folder)):
		split = files.split('.')
		if "log" in split[0]:
			a = tests.single_input(str(config.test_folder) + "/" + files, config.hooks[0])
			b = tests.run_numbers(a)
			status = test.test_run_numbers(default, b, config.run_data, config.run_coef, files)

else:
	
	a = tests.single_input(str(config.test), config.hooks[0])
	b = tests.run_numbers(a)
	status = tests.test_run_numbers(default, b, config.run_data, config.run_coef, config.test)
