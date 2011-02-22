#!/usr/local/bin/python-2.5/bin/python

Info="""
Module name: cfg2lammps.py

Author: (c) Andres Jaramillo-Botero
California Institute of Technology
ajaramil@wag.caltech.edu
Project: pEFF
Version: August 2009

Reads in an eff .cfg file and produces the corresponding lammps data and input files

NOTE: Unsupported functions will be reported in the output log

12/2010: Added support for fixed-core and pseudo-core structures
"""

# import essentials: 
import sys, os 
from math import log10
from shutil import rmtree 
from getopt import gnu_getopt as getopt
import numpy

def printHelp():
    print Info
    print "Usage: python cfg2lammps cfgfile\n"
    return

general="""
# Created %s

# General parameters

variable	sname index %s
log		${sname}.log

units		electron
newton		on
boundary	%s

atom_style      electron

read_data       data.${sname}

pair_style      eff/cut %s
pair_coeff      * *

compute         energies all pair eff/cut
variable        eke equal c_energies[1]
variable        epauli equal c_energies[2]
variable        estatics equal c_energies[3]
variable        errestrain equal c_energies[4]

communicate     single vel yes

compute         peratom all stress/atom
compute         p all reduce sum c_peratom[1] c_peratom[2] c_peratom[3]
variable        press equal -(c_p[1]+c_p[2]+c_p[3])/(3*vol)

compute         effTemp all temp/eff
compute         effPress all pressure effTemp

thermo          %s
thermo_style    custom step etotal pe ke v_eke v_epauli v_estatics v_errestrain temp press v_press 
thermo_modify   temp effTemp press effPress
"""
#%(date,name,boundary,cutoff,period)

minimize="""
# Minimization

min_style       cg
dump            1 %s xyz %s ${sname}.min.xyz
dump            2 %s custom %s ${sname}.min.lammpstrj id type q spin eradius x y z fx fy fz erforce
min_modify      line quadratic
minimize        0 1.0e-5 %s %s

undump		1
undump		2
"""
#%(group,period,group,period,iterations,fcalls)

single_pt="""
# Single point energy

run		0
"""
dynamics="""

# %s Dynamics

timestep        %s

fix             %s

dump            1 %s custom %s ${sname}.%s.lammpstrj id type q spin eradius x y z
dump		2 %s custom %s ${sname}.%s.xyz

run             %s

unfix           1
undump          1
undump		2
"""

task={'single_pt':single_pt,'minimize':minimize,'dynamics':dynamics}

q2m={1:'1.007940',2:'4.002602',3:'6.941000',4:'9.012182',5:'10.811000',6:'12.010700',7:'14.006700',8:'15.999400',
     9:'18.9984032',10:'20.179700',11:'22.98976928',12:'24.305000',13:'26.9815386',14:'28.085500',15:'30.973762',
     16:'32.065000',17:'35.453000',18:'39.948000'}

def generate_lammps_input(infile):

    # Defaults values
    ensemble={"nve":"1 %s nve/eff",'nvt':"1 %s nvt/eff %s %s %s %s",'npt':"1 %s npt/eff %s %s %s %s %s %s"}
    boundary="f f f"
    xbound="-1000.000 1000.0 xlo xhi\n"
    ybound="-1000.000 1000.0 ylo yhi\n"
    zbound="-1000.000 1000.0 zlo zhi\n"
    cutoff=1000.0
    period="1"
    emass=0
    vels=""

    datafile=open("data."+infile[:-4],'w')
    scriptfile=open("in."+infile[:-4],'w')

    print "Reading %s ... [WAIT]"%infile,
    fin = open(infile,'r')
    lines = fin.xreadlines()
    print 7*"\b"+"[DONE]"
    
    numcores=0
    numnuclei=0
    numelec=0
    cores={}
    nuclei={}
    electrons={}
    masses=[]
    massstr="Masses\n\n"
    types=1
    q2type={}
    Tflag=False	# Default ensemble is NVE
    steps='1000'

    print "Extracting run parameters from %s ...  "%(infile),
    for line in lines:

      # 1st level keywords
      if line.find("@params")==0:
        flag='params'
        continue
      elif line.find("@cores")==0:
        flag='cores'
        continue
      elif line.find("@nuclei")==0:
        flag='nuclei'
        continue
      elif line.find("@electrons")==0:
        flag='electrons'
        continue
      elif line.find("@nuc_velocities")==0:
        flag='n_vels'
        continue
      elif line.find("@elec_velocities")==0:
        flag='e_vels'
        continue
      elif line.find("@nuc_masses")==0:
        flag='n_mass'
        continue
      elif line.find("@elec_masses")==0:
        flag='e_mass'
        continue
      elif line.find("@restraints")==0:
        flag='restraints'
        continue

      # 2nd level keywords
      if flag=='params':
        if line.find("calc")>=0:
          op=line.split()[2]
        if line.find("print_every")>=0:
          period=line.split()[2]
        if line.find("num_steps")>=0:
          steps=line.split()[2]
        if line.find("min_freeze")>=0:
          setforce="velocity\t% set 0.0 0.0 0.0\nfix\tfreeze %s setforce 0.0 0.0 0.0"%(line.split()[2],line.split()[2])
        if line.find("thermostat")>=0:
          tflag=True
          #ensemble="fix\t1 all nvt/eff "
        if line.find("start_temperature")>=0:
          Tstart=line.split()[2]
          #ensemble+=Tstart
        if line.find("end_temperature")>=0:
          Tstop=line.split()[2]
          #ensemble+=Tstop
        if line.find("andersen_coupling")>=0 or line.find("nose_hoover_coupling")>=0:
          Tdamp=line.split()[2]
          #ensemble+=Tdamp
        if line.find("dt")>=0:
          dt=line.split()[2]
        if line.find("electron_mass")>=0:
          emass=line.split()[2]
        if line.find("adaptive_step_size")>=0:
          continue
        if line.find("adaptive_energy")>=0:
          continue
        if line.find("e_field_freq")>=0:
          continue        
        if line.find("e_field_packet_duration")>=0:
          continue
        if line.find("e_field")>=0:
          field=line.split()[2:5]
          efield="fix\field all efield %s %s %s"%(field[0],field[1],field[2])
        if line.find("e_field_packet_duration")>=0:
          continue
        if line.find("set_limit")>=0:
          continue	# need to add this contraint
        if line.find("set_limit_stiffness")>=0:
          continue
        if line.find("output_position")>=0:
          dump_pos="dump\t1 all custom %s ${sname}.lammpstrj id type q spin eradius x y z "%(period)
        if line.find("output_velocities")>=0:
          dump_pos+="vx vy vz "
        if line.find("output_energy_forces")>=0:
          dump_pos="compute\tenergy all pe/atom\n"+dump_pos
          dump_pos+="c_energy fx fy fz\n"
        if line.find("output_restart")>=0:
          restart="restart\t%s ${sname}.restart1 ${sname}.restart2"%(period)
        if line.find("output_restraints")>=0:
          continue
        if line.find("ewald_re_cutoff")>=0 or line.find("ewald_autoset")>=0 or line.find("ewald_log_precision")>=0 or line.find("ewald_max_re")>=0 or \
          line.find("ewald_r_cutoff")>=0 or line.find("ewald_k_cutoff")>=0 or line.find("ewald_nuc_r")>=0:
          continue
        if line.find("periodic")>=0: 
          bounds=line.split()[2]
          if bounds=="True": boundary="p p p"
          elif bounds=="minimage_x": boundary="p f f"
          elif bounds=="minimage_xy": boundary="p p f"
          elif bounds=="minimage_y": boundary="f p f"
          elif bounds=="minimage_xyz": boundary="p p p"
          elif bounds=="minimage_z": boundary="f f p"
        if line.find("x_bound")>=0:
          xbnds=line.split()[2:4]
          xbound="%s %s xlo xhi\n"%(xbnds[0],xbnds[1])
        if line.find("y_bound")>=0:
          ybnds=line.split()[2:4]
          ybound="%s %s ylo yhi\n"%(ybnds[0],ybnds[1])
        if line.find("z_bound")>=0:
          zbnds=line.split()[2:4]
          zbound="%s %s zlo zhi\n"%(zbnds[0],zbnds[1])
        if line.find("taper_cutoff")>=0:
          cutoff=line.split()[2]            
        continue

      if flag=='cores' and len(line)>1:
        numcores+=1
        ln=line.split()
        nc=' '.join(ln[0:3])
        q=ln[3]
        spin='3'
        radius=ln[4]
        m=q2m[int(float(q))]
        if m not in masses:
          masses.append(m)
          massstr+="%d %s\n"%(types,m)
          q2type[q]=types
          types+=1
        cores[numcores]=[nc,q,spin,radius]
        continue

      if flag=='nuclei' and len(line)>1:
        numnuclei+=1
        ln=line.split()
        np=' '.join(ln[0:3])
        q=ln[3]
        m=q2m[int(float(q))]
        if m not in masses:
          masses.append(m)
          massstr+="%d %s\n"%(types,m)
          q2type[q]=types
          types+=1
        nuclei[numnuclei]=[np,q]
        continue

      if flag=='electrons' and len(line)>1:
        numelec+=1
        ln=line.split()
        ep=' '.join(ln[0:3])
        spin=ln[3]
        radius=ln[4]
        electrons[numelec]=[ep,spin,radius]
        if numelec==1: 
          if emass!=0: massstr+="%d %s\n\n"%(types,emass)	# electron mass=1
          else: massstr+="%d 1.000000\n\n"%(types)
        continue

      if flag=='n_vels' and len(line)>1:
        vels+=line+" 0.0"
        continue

      if flag=='e_vels' and len(line)>1:
        ln=line.split()
        ln[0]=ln[0]+numnuclei
        vels+=ln[0]+" "+ln[1]+" "+ln[2]+" "+ln[3]+" "+ln[4]+"\n"
        continue

      if flag=='n_mass' and len(line)>1:
        print "Setting nuclear masses is unsupported\n"
        continue

      if flag=='e_mass' and len(line)>1:
        print "Setting electron masses is unsupported\n"
        continue

    print "\bDone"

    # Build data file
    print "Writing datafile to %s ...                   "%('data.'+infile),
    sys.stdout.flush()
    print "\b"*19+"General section  ",
    datafile.writelines("Created using cfg2lammps (c) AJB-2009\n\n%d atoms\n%d atom types\n\n%s%s%s\n"%(numcores+numnuclei+numelec,types,xbound,ybound,zbound))
    print "\b"*19+"Masses section   ",
    datafile.writelines(massstr)
    print "\b"*19+"Atoms section    ",
    datafile.writelines("Atoms\n\n")
    for n in range(numcores):
      datafile.writelines("%d %d %2.2f %s %s %s\n"%(n+1,q2type[cores[n+1][1]],float(cores[n+1][1]),cores[n+1][2],cores[n+1][3],cores[n+1][0]))
    for n in range(numnuclei):
      datafile.writelines("%d %d %2.2f 0 0.0 %s\n"%(n+numcores+1,q2type[nuclei[n+1][1]],float(nuclei[n+1][1]),nuclei[n+1][0]))
    for e in range(numelec):
      datafile.write("%d %d 0.0 %s %s %s\n"%(e+numnuclei+numcores+1,types,electrons[e+1][1],electrons[e+1][2],electrons[e+1][0]))
    print "\b"*19+"Velocities section\n",
    datafile.writelines(vels)
    datafile.writelines("\n")
    print "DONE .... GOODBYE !!"
    datafile.close()

    # Build input script
    import datetime
    scriptfile.writelines(general%(datetime.date.today(),infile[:-4],boundary,cutoff,period))
    if op=='minimize':
      scriptfile.writelines(minimize%('all',period,'all',period,steps,'10000'))
      #%(group,period,group,period,iterations,fcalls)
    elif op=='single_pt':
      scriptfile.writelines(single_pt%())
    elif op=='dynamics': 
      if Tflag==True:
        scriptfile.writelines(dynamics%('NVT',dt,ensemble['nvt']%('all',Tstart,Tstop,Tdamp,''),'all',period,'nvt','all',period,'nve',steps))
        #%(ensemble,dt,group,ensemble%(group,tstart,tstop,tdamp,options))
      else:
        scriptfile.writelines(dynamics%('NVE',dt,ensemble['nve']%('all'),'all',period,'nve','all',period,'nve',steps))
        #%(ensemble,dt,group,ensemble%(group))
    scriptfile.writelines("\n")

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

    generate_lammps_input(infile)
