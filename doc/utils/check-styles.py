#!/usr/bin/env python

from __future__ import print_function
from glob import glob
import os, re

headers = glob(os.path.join('..', '..', 'src', '*', '*.h'))
headers += glob(os.path.join('..', '..', 'src', '*.h'))

angle = {}
atom = {}
body = {}
bond = {}
command = {}
compute = {}
dihedral = {}
dump = {}
fix = {}
improper = {}
integrate = {}
kspace = {}
minimize = {}
pair = {}
reader = {}
region = {}

upper = re.compile("[A-Z]+")
gpu = re.compile("(.+)/gpu$")
intel = re.compile("(.+)/intel$")
kokkos = re.compile("(.+)/kk$")
kokkos_skip = re.compile("(.+)/kk/(host|device)$")
omp = re.compile("(.+)/omp$")
opt = re.compile("(.+)/opt$")
removed = re.compile("(.+)Deprecated$")

def register_style(list,style,info):
    if style in list.keys():
        list[style]['gpu'] += info['gpu']
        list[style]['intel'] += info['intel']
        list[style]['kokkos'] += info['kokkos']
        list[style]['omp'] += info['omp']
        list[style]['opt'] += info['opt']
        list[style]['removed'] += info['removed']
    else:
        list[style] = info

for h in headers:
    # print("Checking ", h)
    fp = open(h)
    text = fp.read()
    matches = re.findall("(.+)Style\((.+),(.+)\)",text,re.MULTILINE)
    for m in matches:

        # skip over internal styles w/o explicit documentation
        style = m[1]
        if upper.match(style):
            continue

        # detect, process, and flag suffix styles:
        info = { 'kokkos':  0, 'gpu':     0, 'intel':   0, \
                 'omp':     0, 'opt':     0, 'removed': 0 }
        suffix = kokkos_skip.match(style)
        if suffix:
            continue
        suffix = gpu.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['gpu'] = 1
        suffix = intel.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['intel'] = 1
        suffix = kokkos.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['kokkos'] = 1
        suffix = omp.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['omp'] = 1
        suffix = opt.match(style)
        if suffix:
            style = suffix.groups()[0]
            info['opt'] = 1
        deprecated = removed.match(m[2])
        if deprecated:
            info['removed'] = 1

        # register style and suffix flags
        if m[0] == 'Angle':
            register_style(angle,style,info)
        elif m[0] == 'Atom':
            register_style(atom,style,info)
        elif m[0] == 'Body':
            register_style(body,style,info)
        elif m[0] == 'Bond':
            register_style(bond,style,info)
        elif m[0] == 'Command':
            register_style(command,style,info)
        elif m[0] == 'Compute':
            register_style(compute,style,info)
        elif m[0] == 'Dihedral':
            register_style(dihedral,style,info)
        elif m[0] == 'Dump':
            register_style(dump,style,info)
        elif m[0] == 'Fix':
            register_style(fix,style,info)
        elif m[0] == 'Improper':
            register_style(improper,style,info)
        elif m[0] == 'Integrate':
            register_style(integrate,style,info)
        elif m[0] == 'KSpace':
            register_style(kspace,style,info)
        elif m[0] == 'Minimize':
            register_style(minimize,style,info)
        elif m[0] == 'Pair':
            register_style(pair,style,info)
        elif m[0] == 'Reader':
            register_style(reader,style,info)
        elif m[0] == 'Region':
            register_style(region,style,info)
        else:
            print("Skipping over: ",m)
fp.close()

