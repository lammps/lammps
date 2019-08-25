# Makefile for the q-TIP4P-F examples
# 
# Copyright (C) 2013, Joshua More and Michele Ceriotti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http.//www.gnu.org/licenses/>.

.PHONY: all clean graphene h2o-piglet_2 h2o-piglet_4 h2o-piglet_8 h2o-rpc h2o-pimd
all: graphene h2o-piglet_2 h2o-piglet_4 h2o-piglet_8 h2o-rpc h2o-pimd

LAMMPS:=../../../../../src/lmp_ubuntu
IPI:=../../../i-pi

define run_lammps
  for i in `seq 1 $1`; do \
    $(LAMMPS) < $2 & \
  done;
endef
-include make.in

graphene:
	cd graphene; $(IPI) input.xml & \
	sleep 20; \
        $(call run_lammps,4,in.kappa.Graphene) \
	wait

h2o-piglet_2:
	cd h2o-piglet.2; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,2,in.water) \
	wait
	
h2o-piglet_4:
	cd h2o-piglet.4; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

h2o-piglet_8:
	cd h2o-piglet.8; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

h2o-pimd:
	cd h2o-pimd; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water) \
	wait

h2o-rpc:
	cd h2o-pimd-rpc; $(IPI) input.xml & sleep 5; \
        $(call run_lammps,4,in.water_shortrange) \
        $(call run_lammps,1,in.water_longrange) \
	wait

clean:
	rm -f */gle_lammps.* */RESTART */EXIT */log.lammps; \
	cd h2o-pimd; \
	rm -f RESTART EXIT log.lammps no_rpc.*; cd ..; \
	cd h2o-pimd-rpc; \
	rm -f RESTART EXIT log.lammps rpc.*; cd ..; \
	cd graphene; \
	rm -f graph.* RESTART EXIT log.lammps graphene.lammpstraj nohup.out; cd ..
