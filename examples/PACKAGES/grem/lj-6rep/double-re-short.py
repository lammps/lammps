#!/usr/bin/env python2.7

import os, sys
from numpy import *
import numpy.random

### Runs replica exchange with gREM (fix grem) for unlimited number of replicas on a set number of processors. This script is inefficient, but necessary if wanting to run with hundreds of replicas on relatively few number of procs.


### read number of processors from the command line
nproc = int(sys.argv[1])

### path to simulation directory
path = os.getcwd()

### path to LAMMPS executable
lmp = sys.argv[2]

### LAMMPS input name
inp = sys.argv[3]

### define pressure for simulations (0 if const V)
pressure = 0

### some constants for gREM, must match with LAMMPS input file!
H = -30000
eta = -0.01
#kB = 0.000086173324 # eV (metal)
kB = 0.0019872 # kcal/mol (real)

### define lambdas - script assumes that there are already existing directories with all files necessary to run
lambdas=[400,405,410,415,420,425]
ll = len(lambdas)

### define number of exchanges
starting_ex = int(loadtxt("lastexchange"))
how_many_ex = 5
max_exchange = starting_ex+how_many_ex

### array with walkers
walker = loadtxt("lastwalker")

### initiate array with enthalpies
enthalpy = zeros(ll)
aver_enthalpy = zeros(ll)

for exchange in arange(starting_ex,max_exchange):
	print "run", exchange
	for l in range(ll):
		#print "replica", l
		os.chdir(path+"/%s" % lambdas[l])
		#os.system("cp restart_file restart_file%d" % exchange)
                if (nproc > 1):
                    os.system("mpirun -np %d " % (nproc) + lmp + " -in ../" + inp + " -var lambda %g -var eta %g -var enthalpy %g > output" % (lambdas[l], eta, H))
                if (nproc == 1):
                    os.system(lmp + " -in ../" + inp + " -var lambda %g -var eta %g -var enthalpy %g > output" % (lambdas[l], eta, H))
		os.system("grep -v '[a-zA-Z]' output | awk '{if(NF==6 && NR>19)print $0}' | awk '{print $3}' >ent")
		enthalpy[l] = os.popen("tail -n 1 ent").read()
		ee = loadtxt("ent")
		aver_enthalpy[l] = mean(ee[-1])
#		os.system("mv dump.dcd dump%d.dcd" % exchange)
		os.system("mv log.lammps log%d.lammps" % exchange)
		os.system("mv final_restart_file final_restart_file%d" % exchange)
		os.system("mv ent ent%d" % exchange)
		os.system("bzip2 log%d.lammps ent%d" % (exchange,exchange))
		os.system("cp final_restart_file%d restart_file" % exchange)

	### replicas will be exchanged based on enthalpy order, not replicas order (termostat order)
	#entalpy_sorted_indices = enthalpy.argsort()
	aver_entalpy_sorted_indices = aver_enthalpy.argsort()

	### choose pair of replicas for exchange attempt based on enthalpy order
	pp = random.random_integers(0,ll-2)
	first = aver_entalpy_sorted_indices[pp]
	second = aver_entalpy_sorted_indices[pp+1]
	#if (first>second):
	#	tmp = first
	#	first = second
	#	second = tmp
	print "pair1:", first, second

	### calculate weights for exchange criterion
	w1 = log(lambdas[first]+eta*(enthalpy[first]-1*H))
	w2 = log(lambdas[first]+eta*(enthalpy[second]-1*H))
	w3 = log(lambdas[second]+eta*(enthalpy[first]-1*H))
	w4 = log(lambdas[second]+eta*(enthalpy[second]-1*H))
	weight = (w4-w3+w1-w2)/eta/kB

	### generate randon number for exchange criterion and calc its log
	LOGRANDNUM = log(random.random())

	### wyzeruj warunki
	compare1 = 0
	compare2 = 0

	if (weight>0):
		compare1 = 1
	if (weight>LOGRANDNUM):
		compare2 = 1

	### exchange restart files if exchange condition is satisfied
	if (compare1>0 or compare2>0):
		print "exchange1 accepted for pair", first, second, lambdas[first], lambdas[second], "with compares as", compare1, compare2, "weight as", weight, "and lograndnum", LOGRANDNUM
		os.system("cp %s/%s/final_restart_file%d %s/%s/restart_file" % (path,lambdas[first],exchange,path,lambdas[second]))
		os.system("cp %s/%s/final_restart_file%d %s/%s/restart_file" % (path,lambdas[second],exchange,path,lambdas[first]))
		### update walkers
		tmp1=walker[first]
		tmp2=walker[second]
		walker[first]=tmp2
		walker[second]=tmp1
	else:
		print "exchange1 not accepted for pair", first, second, lambdas[first], lambdas[second], "with compares as", compare1, compare2, "weight as", weight, "and lograndnum", LOGRANDNUM

	### choose again pair of replicas for exchange attempt based on enthalpy order
	### but make sure this pair is different than the first pair
	if_different = 0
	while if_different<1:
		pp2 = random.random_integers(0,ll-2)
		third = aver_entalpy_sorted_indices[pp2]
		fourth = aver_entalpy_sorted_indices[pp2+1]
		if (third!=first and third!=second and third!=aver_entalpy_sorted_indices[pp-1]):
			if_different = 1

	print "pair2:", third, fourth

	### calculate weights for exchange criterion
	w1 = log(lambdas[third]+eta*(enthalpy[third]-1*H))
	w2 = log(lambdas[third]+eta*(enthalpy[fourth]-1*H))
	w3 = log(lambdas[fourth]+eta*(enthalpy[third]-1*H))
	w4 = log(lambdas[fourth]+eta*(enthalpy[fourth]-1*H))
	weight = (w4-w3+w1-w2)/eta/kB

	### generate randon number for exchange criterion and calc its log
	LOGRANDNUM = log(random.random())

	### wyzeruj warunki
	compare1 = 0
	compare2 = 0

	if (weight>0):
		compare1 = 1
	if (weight>LOGRANDNUM):
		compare2 = 1

	### exchange restart files if exchange condition is satisfied
	if (compare1>0 or compare2>0):
		print "exchange2 accepted for pair", third, fourth, lambdas[third], lambdas[fourth], "with compares as", compare1, compare2, "weight as", weight, "and lograndnum", LOGRANDNUM
		os.system("cp %s/%s/final_restart_file%d %s/%s/restart_file" % (path,lambdas[third],exchange,path,lambdas[fourth]))
		os.system("cp %s/%s/final_restart_file%d %s/%s/restart_file" % (path,lambdas[fourth],exchange,path,lambdas[third]))
		### update walkers
		tmp1=walker[third]
		tmp2=walker[fourth]
		walker[third]=tmp2
		walker[fourth]=tmp1
	else:
		print "exchange2 not accepted for pair", third, fourth, lambdas[third], lambdas[fourth], "with compares as", compare1, compare2, "weight as", weight, "and lograndnum", LOGRANDNUM
	#print "walkers:", walker
	print "".join(["%d " % x for x in walker])
	sys.stdout.flush()

	lastwalker = open(path + "/lastwalker", "w")
	lastwalker.write("".join(["%d " % w for w in walker]))
	lastwalker.close()

	lastexchange = open(path + "/lastexchange", "w")
	lastexchange.write("%d" % (exchange+1))
	lastexchange.close()

