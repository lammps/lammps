#!/usr/local/bin/python-2.5/bin/python

Info="""
Module name: lmp2data.py 

Author: (c) Andres Jaramillo-Botero
California Institute of Technology
ajaramil@wag.caltech.edu
Project: pEFF
Version: August 2009

Extracts the electron radii from a lammps trajectory dump of style custom:

dump    1 all custom period dump_file id type x y z spin radius ...

NOTE: The radius must be the i'th column per trajectory entry in the dump file

"""

# import essentials: 
import sys, os 
from math import log10
from shutil import rmtree 
from getopt import gnu_getopt as getopt
import numpy

def printHelp():
    print Info
    print "Usage: python lmp2data.py test.lammpstrj\n"
    return

def makeradii(infile,outfile,column,flag_all):

    print "Reading %s ... [WAIT]"%infile,
    fin = open(infile,'r')
    lines = fin.xreadlines()
    print 7*"\b"+"[DONE]"
    frame=0
    radii=[]
    # grep the number of frames and atoms/frame
    os.system("grep TIMESTEP %s | wc -l > frames; grep -m 1 -A 1 ATOMS %s > atoms; grep -m 1 \"ITEM: ATOMS\" %s > params"%(infile,infile,infile))
    tmp=open("frames",'r')
    frames=int(tmp.readline().split()[0])
    tmp.close()
    tmp=open("atoms",'r')
    atoms=int(tmp.readlines()[1].split()[0])
    tmp.close()
    tmp=open("params",'r')
    ids=tmp.readline().split()[2:]
    os.system("rm -rf frames atoms params")
    arry=numpy.zeros((atoms,frames),dtype=str)
    framecnt=0
    header=9
    ecount=0
    if flag_all==True: atom_type="nuclei and electron"
    else: atom_type="electron"
    print "Extracting %s %s per frame from %s ...  "%(atom_type,ids[column],infile),
    for i,line in enumerate(lines):
      lo=(atoms+header)*framecnt+header
      hi=lo+atoms
      if (i<lo): 
        continue
      elif (i >= lo) and (i < hi):
        lparse=line.split()
        id=int(lparse[0])
#        r=float(lparse[column-1])
        r=lparse[column]
#        if (float(r)!=0):
        arry[id-1][framecnt]=r
        print arry[id-1][framecnt],r,raw_input()
        if (float(r)!=0) and (framecnt==0): ecount+=1
#        else: arry[id-1][framecnt]=r
        if (i==lo+1): 
          sys.stdout.write("%d/%d%s"%(framecnt+1,frames,(int(log10(framecnt+1))+3+int(log10(frames)))*"\b"))
          sys.stdout.flush()
      if (i == hi+1): 
        framecnt+=1

    print
    if outfile=="": 
      outfile=infile+'.%s'%(ids[column])
      fout=open(outfile,'w')
    else: fout=open(outfile,'w')
    print "Writing %s/frame table to %s ... "%(ids[column],outfile),
    sys.stdout.flush()

    for i in range(frames):
      fout.writelines('\tF'+str(i))
    fout.writelines("\n")
    e=1
    for a in range(atoms):
      if flag_all==True:
        sys.stdout.write("%d/%d%s"%(a+1,atoms,(int(log10(a+1))+int(log10(atoms))+3)*"\b"))
        sys.stdout.flush()
        fout.writelines("%d\t"%(a+1))
        for f in range(frames):
          fout.writelines("%s\t"%(arry[a][f]))
        fout.writelines("\n")
      else:
        if arry[a][0] == 0.0:
          continue
        else:
          sys.stdout.write("%d/%d%s"%(e,ecount,(int(log10(e))+int(log10(ecount))+3)*"\b"))
          sys.stdout.flush()
          e+=1
          fout.writelines("%d\t"%(a+1))
          for f in range(frames):
            fout.writelines("%s\t"%(arry[a][f]))
          fout.writelines("\n")

    print
    print "DONE .... GOODBYE !!"
    fout.close()
    fin.close()

if __name__ == '__main__':

    # set defaults
    outfile = ""
    flag_all = False
    column=6	# default = radius

    # check for input:
    opts, argv = getopt(sys.argv[1:], 'c:o:ha')

    # if no input, print help and exit
    if len(argv) != 1:
        printHelp()
        sys.exit(1)
    else:
        infile=argv[0]

    # read options
    for opt, arg in opts:
        if opt == '-h':             # -h: print help
          printHelp()
        if opt == '-o':		# output file name
          outfile=arg
        if opt == '-a':		# all nuclii+electrons
          flag_all=True
        if opt == '-c':		# select column from lammpstrj file to tabulate
          column=int(arg)

    makeradii(infile,outfile,column,flag_all)
