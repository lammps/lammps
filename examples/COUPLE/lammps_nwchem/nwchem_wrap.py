#!/usr/bin/env python

# ----------------------------------------------------------------------
# LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
# https://www.lammps.org/ Sandia National Laboratories
# Steve Plimpton, sjplimp@sandia.gov
# ----------------------------------------------------------------------

# Syntax: nwchem_wrap.py file/zmq ao/pw input_template
# file/zmg = messaging mode, must match LAMMPS messaging mode
# ao/pw = basis set mode, selects between atom-centered and plane-wave
#         the input_template file must correspond to the appropriate basis set mode:
#         the "ao" mode supports the scf and dft modules in NWChem,
#         the "pw" mode supports the nwpw module.
# input_template = NWChem input file used as template, must include a
#                  "geometry" block with the atoms in the simulation, dummy
#                  xyz coordinates should be included (but are not used).
#                  Atom ordering must match LAMMPS input.

# wrapper on NWChem
#   receives message with list of coords
#   creates NWChem inputs
#   invokes NWChem to calculate self-consistent energy of that config
#   reads NWChem outputs
#   sends message with energy, forces, pressure to client

from __future__ import print_function
import sys

version = sys.version_info[0]
if version == 3:
  sys.exit("The CSlib python wrapper does not yet support python 3")

import subprocess
import re
import os
import shutil
from cslib import CSlib

# comment out 2nd line once 1st line is correct for your system

nwchemcmd = "mpirun -np 1 /usr/bin/nwchem"
nwchemcmd = "touch tmp"

# enums matching FixClientMD class in LAMMPS

SETUP,STEP = range(1,2+1)
DIM,PERIODICITY,ORIGIN,BOX,NATOMS,NTYPES,TYPES,COORDS,UNITS,CHARGE = range(1,10+1)
FORCES,ENERGY,VIRIAL,ERROR = range(1,4+1)

# -------------------------------------
# functions

# error message and exit

def error(txt):
  print("ERROR:",txt)
  sys.exit(1)

# -------------------------------------
# read initial input file to setup problem
# return natoms

def nwchem_setup_ao(input):

  template = open(input,'r')

  geometry_block = False
  natoms = 0

  while True:
    line = template.readline()
    if not line: break

    if geometry_block and re.search("end",line):
      geometry_block = False
    if geometry_block and not re.match("#",line) :
      natoms += 1
    if re.search("geometry",line):
      geometry_block = True

  return natoms

# -------------------------------------
# write a new input file for NWChem
# assumes the NWChem input geometry is to be specified in angstroms

def nwchem_input_write_ao(input,coords):

  template = open(input,'r')
  new_input = open("nwchem_lammps.nw",'w')

  geometry_block = False
  i = 0

  while True:
    line = template.readline()
    if not line: break

    if geometry_block and not re.match("#",line) and re.search("end",line):
      geometry_block = False
      if os.path.exists("nwchem_lammps.movecs"):
        # The below is hacky, but one of these lines will be ignored
        # by NWChem depending on if the input file is for scf/dft.
        append = "\nscf\n vectors input nwchem_lammps.movecs\nend\n"
        append2 = "\ndft\n vectors input nwchem_lammps.movecs\nend\n"
        line = line + append + append2

    if geometry_block and not re.match("#",line):
      x = coords[3*i+0]
      y = coords[3*i+1]
      z = coords[3*i+2]
      coord_string = "  %g %g %g \n" % (x,y,z)
      atom_string = line.split()[0]
      line = atom_string + coord_string
      i += 1

    if (not re.match("#",line)) and re.search("geometry",line):
      geometry_block = True
      line = "geometry units angstrom noautosym nocenter\n"

    print(line,file=new_input,end='')

  new_input.close()

# -------------------------------------
# read a NWChem output nwchem_lammps.out file

def nwchem_read_ao(natoms, log):

  nwchem_output = open(log, 'r')
  energy_pattern = r"Total \w+ energy"
  gradient_pattern = "x          y          z           x          y          z"

  eout = 0.0
  fout = []

  while True:
    line = nwchem_output.readline()
    if not line: break

    # pattern match for energy
    if re.search(energy_pattern,line):
      eout = float(line.split()[4])

    # pattern match for forces
    if re.search(gradient_pattern, line):
      for i in range(natoms):
        line = nwchem_output.readline()
        forces = line.split()
        fout += [float(forces[5]), float(forces[6]), float(forces[7])]

  # convert units
  hartree2eV = 27.21138602
  bohr2angstrom = 0.52917721092
  eout = eout * hartree2eV
  fout = [i * -hartree2eV/bohr2angstrom for i in fout]

  return eout,fout

# -------------------------------------
# read initial planewave input file to setup problem
# return natoms,box
def nwchem_setup_pw(input):

  template = open(input,'r')

  geometry_block = False
  system_block = False
  coord_pattern = r"^\s*\w{1,2}(?:\s+-?(?:\d+.?\d*|\d*.?\d+)){3}"
  natoms = 0
  box = []

  while True:
    line = template.readline()
    if not line: break

    if geometry_block and re.search("system crystal",line):
      system_block = True
      for i in range(3):
        line = template.readline()
        line = re.sub(r'd|D', 'e', line)
        box += [float(line.split()[1])]

    if geometry_block and not system_block and re.match("#",line) and re.search("end",line):
      geometry_block = False

    if system_block and re.search("end",line):
      system_block = False

    if geometry_block and not re.match("#",line) and re.search(coord_pattern,line):
      natoms += 1

    if re.search("geometry",line) and not re.match("#",line):
      geometry_block = True

  return natoms,box

# -------------------------------------
# write a new planewave input file for NWChem
# assumes the NWChem input geometry is to be specified fractional coordinates

def nwchem_input_write_pw(input,coords,box):

  template = open(input,'r')
  new_input = open("nwchem_lammps.nw",'w')

  writing_atoms = False
  geometry_block = False
  system_block = False
  coord_pattern = r"^\s*\w{1,2}(?:\s+-?(?:\d+.?\d*|\d*.?\d+)){3}"
  i = 0

  while True:
    line = template.readline()
    if not line: break

    if geometry_block and re.search("system crystal",line):
      system_block = True

    if geometry_block and not system_block and not re.match("#",line) and re.search("end",line):
      geometry_block = False
      if os.path.exists("nwchem_lammps.movecs"):
        append = "\nnwpw\n vectors input nwchem_lammps.movecs\nend\n"
        line = line + append

    if system_block and re.search("end",line):
      system_block = False

    if geometry_block and not re.match("#",line) and re.search(coord_pattern,line):
      x = coords[3*i+0] / box[0]
      y = coords[3*i+1] / box[1]
      z = coords[3*i+2] / box[2]
      coord_string = "  %g %g %g \n" % (x,y,z)
      atom_string = line.split()[0]
      line = atom_string + coord_string
      i += 1

    if re.search("geometry",line) and not re.match("#",line):
      geometry_block = True

    print(line,file=new_input,end='')

  new_input.close()

# -------------------------------------
# read a NWChem output nwchem_lammps.out file for planewave calculation

def nwchem_read_pw(log):
  nw_output = open(log, 'r')

  eout = 0.0
  sout = []
  fout = []
  reading_forces = False

  while True:
    line = nw_output.readline()
    if not line: break

    # pattern match for energy
    if re.search("PSPW energy",line):
      eout = float(line.split()[4])

    # pattern match for forces
    if re.search("C\.O\.M", line):
      reading_forces = False
    if reading_forces:
      forces = line.split()
      fout += [float(forces[3]), float(forces[4]), float(forces[5])]
    if re.search("Ion Forces",line):
      reading_forces = True

    # pattern match for stress
    if re.search("=== total gradient ===",line):
      stensor = []
      for i in range(3):
        line = nw_output.readline()
        line = line.replace("S =","   ")
        stress = line.split()
        stensor += [float(stress[1]), float(stress[2]), float(stress[3])]
      sxx = stensor[0]
      syy = stensor[4]
      szz = stensor[8]
      sxy = 0.5 * (float(stensor[1]) + float(stensor[3]))
      sxz = 0.5 * (stensor[2] + stensor[6])
      syz = 0.5 * (stensor[5] + stensor[7])
      sout = [sxx,syy,szz,sxy,sxz,syz]

  # convert units
  hartree2eV = 27.21138602
  bohr2angstrom = 0.52917721092
  austress2bar = 294210156.97
  eout = eout * hartree2eV
  fout = [i * hartree2eV/bohr2angstrom for i in fout]
  sout = [i * austress2bar for i in sout]

  return eout,fout,sout

# -------------------------------------
# main program

# command-line args
#
if len(sys.argv) != 4:
  print("Syntax: python nwchem_wrap.py file/zmq ao/pw input_template")
  sys.exit(1)

comm_mode = sys.argv[1]
basis_type = sys.argv[2]
input_template = sys.argv[3]

if comm_mode == "file": cs = CSlib(1,comm_mode,"tmp.couple",None)
elif comm_mode == "zmq": cs = CSlib(1,comm_mode,"*:5555",None)
else:
  print("Syntax: python nwchem_wrap.py file/zmq")
  sys.exit(1)


natoms = 0
box = []
if basis_type == "ao":
  natoms = nwchem_setup_ao(input_template)
elif basis_type == "pw":
  natoms,box = nwchem_setup_pw(input_template)

# initial message for AIMD protocol

msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
if msgID != 0: error("Bad initial client/server handshake")
protocol = cs.unpack_string(1)
if protocol != "md": error("Mismatch in client/server protocol")
cs.send(0,0)

# endless server loop

i = 0
if not os.path.exists("nwchem_logs"):
  os.mkdir("nwchem_logs")

while 1:

  # recv message from client
  # msgID = 0 = all-done message

  msgID,nfield,fieldID,fieldtype,fieldlen = cs.recv()
  if msgID < 0: break

  # SETUP receive at beginning of each run
  # required fields: DIM, PERIODICITY, ORIGIN, BOX,
  #                  NATOMS, COORDS
  # optional fields: others in enum above, but NWChem ignores them

  if msgID == SETUP:

    origin = []
    box_lmp = []
    natoms_recv = ntypes_recv = 0
    types = []
    coords = []

    for field in fieldID:
      if field == DIM:
        dim = cs.unpack_int(DIM)
        if dim != 3: error("NWChem only performs 3d simulations")
      elif field == PERIODICITY:
        periodicity = cs.unpack(PERIODICITY,1)
        if basis_type == "ao":
          if periodicity[0] or periodicity[1] or periodicity[2]:
            error("NWChem AO basis wrapper only currently supports fully aperiodic systems")
        elif basis_type == "pw":
          if not periodicity[0] or not periodicity[1] or not periodicity[2]:
            error("NWChem PW basis wrapper only currently supports fully periodic systems")
      elif field == ORIGIN:
        origin = cs.unpack(ORIGIN,1)
      elif field == BOX:
        box_lmp = cs.unpack(BOX,1)
        if (basis_type == "pw"):
          if (box[0] != box_lmp[0] or box[1] != box_lmp[4] or box[2] != box_lmp[8]):
            error("NWChem wrapper mismatch in box dimensions")
      elif field == NATOMS:
        natoms_recv = cs.unpack_int(NATOMS)
        if natoms != natoms_recv:
          error("NWChem wrapper mismatch in number of atoms")
      elif field == COORDS:
        coords = cs.unpack(COORDS,1)

    if not origin or not box_lmp or not natoms or not coords:
      error("Required NWChem wrapper setup field not received");

  # STEP receive at each timestep of run or minimization
  # required fields: COORDS
  # optional fields: ORIGIN, BOX

  elif msgID == STEP:

    coords = []

    for field in fieldID:
      if field == COORDS:
        coords = cs.unpack(COORDS,1)

    if not coords: error("Required NWChem wrapper step field not received");

  else: error("NWChem wrapper received unrecognized message")

  # unpack coords from client
  # create NWChem input

  if basis_type == "ao":
    nwchem_input_write_ao(input_template,coords)
  elif basis_type == "pw":
    nwchem_input_write_pw(input_template,coords,box)

  # invoke NWChem

  i += 1
  log = "nwchem_lammps.out"
  archive = "nwchem_logs/nwchem_lammps" + str(i) + ".out"
  cmd = nwchemcmd + " nwchem_lammps.nw > " + log
  print("\nLaunching NWChem ...")
  print(cmd)
  subprocess.check_output(cmd,stderr=subprocess.STDOUT,shell=True)

  shutil.copyfile(log,archive)

  # process NWChem output

  if basis_type == "ao":
    energy,forces = nwchem_read_ao(natoms,log)
    virial = [0,0,0,0,0,0]
  elif basis_type == "pw":
    energy,forces,virial = nwchem_read_pw(log)

  # return forces, energy to client
  cs.send(msgID,3)
  cs.pack(FORCES,4,3*natoms,forces)
  cs.pack_double(ENERGY,energy)
  cs.pack(VIRIAL,4,6,virial)

# final reply to client

cs.send(0,0)

# clean-up

del cs
