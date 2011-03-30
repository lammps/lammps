#!/usr/local/bin/python-2.5/bin/python

Info="""
Module name: lmp2radii.py 

Author: (c) Andres Jaramillo-Botero
California Institute of Technology
ajaramil@wag.caltech.edu
Project: pEFF
Version: August 2009

Extracts the electron radii from a lammps trajectory dump of style custom:

dump    1 all custom period dump_file id type q spin eradius x y z...

NOTE: The radius must be the 5th column per trajectory entry in the dump file

"""

# import essentials: 
import sys, os 
from math import log10
from shutil import rmtree 
from getopt import gnu_getopt as getopt
import numpy

def printHelp():
    print Info
    print "Usage: python lmp2radii.pyx test.lammpstrj\n"
    return

def makeradii(infile):

    print "Reading %s ... [WAIT]"%infile,
    fin = open(infile,'r')
    lines = fin.xreadlines()
    print 7*"\b"+"[DONE]"
    frame=0
    radii=[]
    # grep the number of frames and atoms/frame
    os.system("grep TIMESTEP %s | wc -l > frames; grep -m 1 -A 1 ATOMS %s > atoms"%(infile,infile))
    tmp=open("frames",'r')
    frames=int(tmp.readline().split()[0])
    tmp.close()
    tmp=open("atoms",'r')
    atoms=int(tmp.readlines()[1].split()[0])
    tmp.close()
    os.system("rm -rf frames atoms lines")
    arry=numpy.zeros((atoms,frames),dtype=float)
    framecnt=0
    header=9
    ecount=0
    print "Extracting electron radii per frame from %s ...  "%(infile),
    for i,line in enumerate(lines):
      lo=(atoms+header)*framecnt+header
      hi=lo+atoms
      if (i<lo): 
        continue
      elif (i >= lo) and (i < hi):
        lparse=line.split()
        id=int(lparse[0])
        r=float(lparse[4])
        if (r!=0.0): 
          arry[id-1][framecnt]=r
          if (framecnt==0): ecount+=1
        if (i==lo+1): 
          sys.stdout.write("%d/%d%s"%(framecnt+1,frames,(int(log10(framecnt+1))+3+int(log10(frames)))*"\b"))
          sys.stdout.flush()
      if (i == hi+1): 
        framecnt+=1
    print
    print "Writing radii/frame table to %s ... "%(infile+'.out'),
    sys.stdout.flush()
    fout=open(infile+'.out','w')
    for i in range(frames):
      fout.writelines('\tF'+str(i))
    fout.writelines("\n")
    e=1
    for a in range(atoms):
      if arry[a][0] == 0.0: continue
      else:
        sys.stdout.write("%d/%d%s"%(e,ecount,(int(log10(e))+int(log10(ecount))+3)*"\b"))
        sys.stdout.flush()
        e+=1
        fout.writelines("%d\t"%(a+1))
        for f in range(frames):
          fout.writelines("%f\t"%(arry[a][f]))
        fout.writelines("\n")
    print
    print "Done !! (generated radii/frame table) \n"
    fout.close()
    fin.close()

if __name__ == '__main__':

    # set defaults
    
    # check for input:
    opts, argv = getopt(sys.argv[1:], 'h')

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

    makeradii(infile)
