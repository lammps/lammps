#!/usr/bin/env python3
"""
Find all example input files containing commands changed in this branch versus develop.
Companion script to run_tests.py regression tester.
"""

import os, re, sys
from glob import glob
import subprocess

# infer top level lammps dir
LAMMPS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

# get list of changed files relative to the develop branch from git
output = None
try:
    output = subprocess.run('git diff --diff-filter=MA --name-status develop',
                            shell=True, capture_output=True)
except:
    pass

# collect header files to check for styles
headers = []
if output:
    for changed in output.stdout.decode().split():
        if (changed == 'A') or (changed == 'M'): continue
        if not changed.startswith('src/'): continue
        if changed.endswith('.h'): headers.append(changed)
        if changed.endswith('.cpp'): headers.append(changed.replace('.cpp','.h'))

# now loop over header files, search for XxxxStyle() macros and append style name
command = []
atom = []
compute = []
fix = []
pair = []
body = []
bond = []
angle = []
dihedral = []
improper = []
kspace = []
dump = []
region = []
integrate = []
minimize = []

style_pattern = re.compile(r"(.+)Style\((.+),(.+)\)")
upper = re.compile("[A-Z]+")
gpu = re.compile("(.+)/gpu$")
intel = re.compile("(.+)/intel$")
kokkos = re.compile("(.+)/kk$")
kokkos_skip = re.compile("(.+)/kk/(host|device)$")
omp = re.compile("(.+)/omp$")
opt = re.compile("(.+)/opt$")
removed = re.compile("(.*)Deprecated$")

for file in headers:
    with open(file) as f:
        for line in f:
            matches = style_pattern.findall(line)
            for m in matches:
                # skip over internal styles w/o explicit documentation
                style = m[1]
                if upper.match(style):
                    continue

                # skip over suffix styles:
                suffix = kokkos_skip.match(style)
                if suffix:
                    continue
                suffix = gpu.match(style)
                if suffix:
                    continue
                suffix = intel.match(style)
                if suffix:
                    continue
                suffix = kokkos.match(style)
                if suffix:
                    continue
                suffix = omp.match(style)
                if suffix:
                    continue
                suffix = opt.match(style)
                if suffix:
                    continue
                deprecated = removed.match(m[2])
                if deprecated:
                    continue

                # register style and suffix flags
                if m[0] == 'Angle':
                    angle.append(style)
                elif m[0] == 'Atom':
                    atom.append(style)
                elif m[0] == 'Body':
                    register_style(body,style,info)
                elif m[0] == 'Bond':
                    bond.applend(style)
                elif m[0] == 'Command':
                    command.append(style)
                elif m[0] == 'Compute':
                    compute.append(style)
                elif m[0] == 'Dihedral':
                    dihedral.append(style)
                elif m[0] == 'Dump':
                    dump.append(style)
                elif m[0] == 'Fix':
                    fix.append(style)
                elif m[0] == 'Improper':
                    improper.append(style)
                elif m[0] == 'Integrate':
                    integrate.append(style)
                elif m[0] == 'KSpace':
                    kspace.append(style)
                elif m[0] == 'Minimize':
                    minimize.append(style)
                elif m[0] == 'Pair':
                    pair.append(style)
                elif m[0] == 'Region':
                    region.append(style)
                else:
                    pass

if len(command):
    print("Commands: ", '|'.join(command))
if len(atom):
    print("Atom styles: ", '|'.join(atom))
if len(compute):
    print("Compute styles: ", '|'.join(compute))
if len(fix):
    print("Fix styles: ", '|'.join(fix))
if len(pair):
    print("Pair styles: ", '|'.join(pair))
if len(body):
    print("Body styles: ", '|'.join(body))
if len(bond):
    print("Bond styles: ", '|'.join(bond))
if len(angle):
    print("Angle styles: ", '|'.join(angle))
if len(dihedral):
    print("Dihedral styles: ", '|'.join(dihedral))
if len(improper):
    print("Improper styles: ", '|'.join(improper))
if len(kspace):
    print("Kspace styles: ", '|'.join(kspace))
if len(dump):
    print("Dump styles: ", '|'.join(dump))
if len(region):
    print("Region styles: ", '|'.join(region))
if len(integrate):
    print("Integrate styles: ", '|'.join(integrate))
if len(minimize):
    print("Minimize styles: ", '|'.join(minimize))
