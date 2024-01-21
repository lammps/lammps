#!/usr/bin/env python3

from yaml import load
import subprocess
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

runs = subprocess.check_output('gh api repos/lammps/lammps/actions/runs',shell=True)
data = load(runs,Loader=Loader)
while data['total_count'] > 3:
    print('remaining: ', data['total_count'])
    num=1
    for d in data['workflow_runs']:
        print(num, d['id'],d['name'],d['run_number'])
        num += 1
        if num > 4:
            subprocess.call('gh api -X DELETE repos/lammps/lammps/actions/runs/' + str(d['id']), shell=True)
            #print('gh api -X DELETE repos/lammps/lammps/actions/runs/' + str(d['id']))
        else:
            print('skip')
    runs = subprocess.check_output('gh api repos/lammps/lammps/actions/runs',shell=True)
    data = load(runs,Loader=Loader)
