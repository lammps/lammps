# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# vmd tool

# Minimalistic VMD embedding for Pizza.py
# (c) 2010 Axel Kohlmeyer <akohlmey@gmail.com>
# This class will replace the VMD startup script,
#   open a pipe to the executable,
#   and feed it Tcl command lines one at a time

oneline = "Control VMD from python"

docstr = """
v = vmd()		       start up VMD
v.stop()		       shut down VMD instance
v.clear()		       delete all visualizations

v.rep(style)		       set default representation style. One of
			       (Lines|VDW|Licorice|DynamicBonds|Points|CPK) 
v.new(file[,type])   	       load new file (default file type 'lammpstrj')
v.data(file[,atomstyle])       load new data file (default atom style 'full')
v.replace(file[,type])	       replace current frames with new file
v.append(file[,type]) 	       append file to current frame(s)
v.set(snap,x,y,z,(True|False)) set coordinates from a pizza.py snapshot to new or current frame

v.frame(frame)		       set current frame
v.flush()		       flush pending input to VMD and update GUI
v.read(file)		       read Tcl script file (e.g. saved state)
    
v.enter()		       enter interactive shell
v.debug([True|False])	       display generated VMD script commands?
"""

# History
#   11/10, Axel Kohlmeyer (Temple U): original version

# Imports and external programs

import types, os
import numpy

try: from DEFAULTS import PIZZA_VMDNAME
except: PIZZA_VMDNAME = "vmd"
try: from DEFAULTS import PIZZA_VMDDIR
except: PIZZA_VMDDIR = "/usr/local/lib/vmd"
try: from DEFAULTS import PIZZA_VMDDEV
except: PIZZA_VMDDEV = "win"
try: from DEFAULTS import PIZZA_VMDARCH
except: PIZZA_VMDARCH = "LINUX"

# try these settings for a Mac
#PIZZA_VMDNAME = "vmd"
#PIZZA_VMDDIR = "/Applications/VMD\ 1.8.7.app/Contents/vmd"
#PIZZA_VMDDEV = "win"
#PIZZA_VMDARCH = "MACOSXX86"

try: import pexpect
except: 
  print "pexpect from http://pypi.python.org/pypi/pexpect", \
      "is required for vmd tool"
  raise

# Class definition

class vmd:
  
  # --------------------------------------------------------------------

  def __init__(self):
    self.vmddir = PIZZA_VMDDIR
    self.vmdexe = PIZZA_VMDDIR + '/' + PIZZA_VMDNAME + '_' + PIZZA_VMDARCH
    # these are all defaults copied from the vmd launch script
    os.environ['VMDDIR'] = PIZZA_VMDDIR
    os.environ['VMDDISPLAYDEVICE'] = PIZZA_VMDDEV
    os.environ['VMDSCRPOS']    = "596 190"
    os.environ['VMDSCRSIZE']   = "669 834"
    os.environ['VMDSCRHEIGHT'] = "6.0"
    os.environ['VMDSCRDIST']   = "-2.0"
    os.environ['VMDTITLE']     = "on"
    os.environ['TCL_LIBRARY'] = PIZZA_VMDDIR + "/scripts/tcl"
    os.environ['STRIDE_BIN']  = PIZZA_VMDDIR + "/stride_" + PIZZA_VMDARCH
    os.environ['SURF_BIN']    = PIZZA_VMDDIR + "/surf_" + PIZZA_VMDARCH
    os.environ['TACHYON_BIN'] = PIZZA_VMDDIR + "/tachyon_" + PIZZA_VMDARCH
    ldpath = os.environ.get('LD_LIBRARY_PATH','')
    if ldpath == '':
      os.environ['LD_LIBRARY_PATH'] = PIZZA_VMDDIR
    else:
      os.environ['LD_LIBRARY_PATH'] = ldpath + ':' + PIZZA_VMDDIR
    ldpath = os.environ.get('LD_LIBRARY_PATH','')
    if ldpath == '':
      os.environ['PYTHONPATH'] = PIZZA_VMDDIR
    else:
      os.environ['PYTHONPATH'] = PIZZA_VMDDIR + "/scripts/python"
    self.debugme = False
    # open pipe to vmd and wait until we have a prompt
    self.VMD = pexpect.spawn(self.vmdexe)
    self.VMD.expect('vmd >')
              
  # --------------------------------------------------------------------
  # post command to vmd and wait until the prompt returns.
  def __call__(self,command):
    if self.VMD.isalive():
      self.VMD.sendline(command)
      self.VMD.expect('vmd >')
      if self.debugme:
        print "call+result:"+self.VMD.before
    return
    
  # --------------------------------------------------------------------
  # exit VMD
  def stop(self):
    self.__call__("quit")
    del self.VMD

  # --------------------------------------------------------------------
  # force VMD display and GUI update.
  def flush(self):
    self.__call__('display update ui')

  # --------------------------------------------------------------------
  # turn on debugging info
  def debug(self,status=True):
    if status and not self.debugme:
      print 'Turning vmd.py debugging ON.'
    if not status and self.debugme:
      print 'Turning vmd.py debugging OFF.'
    self.debugme = status

  # --------------------------------------------------------------------
  # emulate a regular tcl command prompt
  def enter(self,mode='tcl'):
    self.__call__('menu main off')
    self.__call__('menu main on')
    while 1:
      try:
        command = raw_input("vmd > ")
      except EOFError:
        print "(EOF)"
        self.__call__('menu main off')
        return
      if command == "quit" or command == "exit":
        self.__call__('menu main off')
        return
      if command == "gopython":
        print "gopython not supported here"
        continue
      self.__call__(command)

  # --------------------------------------------------------------------
  # read and execute tcl script file (e.g. a saved state)
  def read(self,filename):
    self.__call__('play ' + filename)
    self.flush()

  # --------------------------------------------------------------------
  # remove all molecules, data and visualizations
  def clear(self):
    self.__call__("mol delete all")

  # --------------------------------------------------------------------
  # navigate to a given frame
  def rep(self,style='Lines'):
    if style == 'Lines' or style == 'VDW' or style == 'Licorice' \
       or style == 'DynamicBonds' or style == 'Points' or style == 'CPK':
      self.__call__('mol default style ' + style)
  # --------------------------------------------------------------------
  # navigate to a given frame
  def frame(self,framespec):
    self.__call__('animate goto ' + str(framespec))

  # --------------------------------------------------------------------
  # load a new molecule from a file supported by a molfile plugin
  def new(self,filename,filetype='lammpstrj'):
    self.__call__('mol new ' + filename + ' type ' + filetype + ' waitfor all')
    self.flush()

  # --------------------------------------------------------------------
  # load a new molecule from a data file via the topotools plugin
  def data(self,filename,atomstyle='full'):
    self.__call__('package require topotools 1.0')
    self.__call__('topo readlammpsdata ' + filename + ' ' + atomstyle)
    self.flush()

  # --------------------------------------------------------------------
  # append all frames from a given file to the current molecule
  def append(self,filename,filetype='lammpstrj'):
    self.__call__('set tmol [molinfo top]')
    self.__call__('array set viewpoints {}')
    self.__call__('foreach mol [molinfo list] { set viewpoints($mol) [molinfo $mol get { center_matrix rotate_matrix scale_matrix global_matrix}]}')
    self.__call__('mol addfile ' + filename + ' mol $tmol type ' + filetype + ' waitfor all')
    self.__call__('foreach mol [molinfo list] { molinfo $mol set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($mol)}')
    self.flush()
    
  # --------------------------------------------------------------------
  # replace all frames of a molecule with those from a given file
  def update(self,filename,filetype='lammpstrj'):
    self.__call__('set tmol [molinfo top]')
    self.__call__('array set viewpoints {}')
    self.__call__('foreach mol [molinfo list] {set viewpoints($mol) [molinfo $mol get { center_matrix rotate_matrix scale_matrix global_matrix}]}')
    self.__call__('animate delete all $tmol')
    self.__call__('mol addfile ' + filename + ' mol $tmol type ' + filetype + ' waitfor all')
    self.__call__('foreach mol [molinfo list] {molinfo $mol set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($mol)}')
    self.flush()
    
  # --------------------------------------------------------------------
  # add or overwrite coordinates with coordinates in a snapshot
  def set(self,snap,x,y,z,append=True):
    self.__call__('set vmdsel [atomselect top all]')
    if append:
      self.__call__('animate dup [molinfo top]')
    cmd = '$vmdsel set {x y z} {'
    for idx in range(0,snap.natoms):
      cmd += ' {'+str(snap[idx,x])+' '+str(snap[idx,y])+' '+str(snap[idx,z])+'}'
    cmd += '}'
    self.__call__(cmd)
    self.__call__('$vmdsel delete ; unset vmdsel')
    self.flush()
