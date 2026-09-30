# QEq cannot be run together with minimize, because each new calculation of charges alters the energy
# hypersurface on which we are looking for a minimum. Minimization frequently gets stuck with the
# message "linesearch alpha is zero" before reaching zero pressure.
#
# Below is a self-consistency loop that minimizes the potential energy while distributing atomic
# charges, which replaces the minimize command. It works as follows:
#   1. For given atomic positions, run QEq to compute charges. The initial charges are q0 and the
#      final ones q.
#   2. Fix charges (q) and run energy minimization with changing box shape to get new atomic
#      coordinates. This gives zero pressure.
#   3. If max|q-q0| >= 0.01, set q0 := q and go to 1.
#      If max|q-q0| < 0.01, we have equilibrium atomic positions and charges => exit.
#
# Local variables end with _ to avoid interference with the main script.

variable dq_ atom abs(q-f_q0_)
compute dqmax_ all reduce max v_dq_
variable dqmax_ equal c_dqmax_

variable iter loop 100   # number of iterations to equilibrate charges & get zero pressure
label self_consistency
  fix q0_ all store/state 0 q
  fix fq all qeq/slater 1 12.0 1e-6 1000 coul/streitz alpha 0.3 dsf 2.0
  run 0	      	      	      	# calculate charges -> nonzero pressure

  unfix fq
  minimize 0 1e-6 10000 10000 	# fix charges -> zero pressure

  print "***** self-consistent loop: iter=${iter}, dqmax=${dqmax_} *****"
  if "${dqmax_} < 0.01" then "jump SELF break"
  next iter
jump SELF self_consistency

label break
unfix q0_
uncompute dqmax_
variable iter delete
