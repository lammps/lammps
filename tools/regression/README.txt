regression.py

Tool for numerical comparisions of benchmark log files
Created by Stan Moore (SNL), email: stamoor at sandia.gov
based on benchmark.py created by Reese Jones (SNL)
Requires log.py from Pizza.py toolkit, http://pizza.sandia.gov/


SYNTAX: regression.py <descriptor> <LAMMPS_args> <test_dirs> <options>

    descriptor = any string without spaces, appended to log files, allows multiple tests in the same directory
    LAMMPS_args = string to launch the benchmark calculation
      the path to the executable must be an absolute path
      e.g. ~/lammps/src/lmp_g++ or "mpirun -np 4 ~/lammps/src/lmp_g++ -v x 10"
    test_dirs = list of one or more dirs to recursively search for scripts
      scripts = any in.* file
    options = one or more keyword/value pairs
      wildcards are expanded by Python, not the shell, so it may be necessary to escape * as \*

"Gold standard" logfiles are automatically generated if they don't exist


EXAMPLES:

python regression.py mpi_16 "mpiexec -np 16 ~/lammps/src/lmp_mpi" ~/regression/examples -only colloide 2>&1 |tee test_mpi_16.out
python regression.py kk_gpu_2 "mpiexec -np 2 .lmp_kokkos_cuda_openmpi -k on g 2 -sf kk -pk kokkos comm/forward device comm/exchange device newton on neigh half" ~/lammps/examples -error_norm L1 -relative_error True -tolerance 0.05 -min-same-rows 2 -exclude hugoniostat nb3b colloid indent snap peri dreiding ASPHERE pour streitz srd min balance USER/cg-cmm/sds-monolayer USER/gauss-diel USER/awpmd/H USER/fep USER/misc/basal USER/lb USER/eff USER/drude USER/qtb/alpha_quartz_qbmsst USER/tally ELASTIC tests/from_Jeff ellipse body dipole comb kim tests/from_Dan rigid VISCOSITY micelle USER/gle USER/qtb/methane_qbmsst USER/qtb/methane_qtb USER/cg-cmm/peg-verlet USER/pimd/para-h2 reax voronoi USER/qtb/alpha_quartz_qtb MC meam shear tests/GCMC USER/dpd 2>&1 |tee test_kk_gpu_2_half.out


OPTIONS:
      -exclude <subdir1 subdir2* ...>
        do not run tests from these sub-dirs or their children
        default = none
      -only <subdir1 subdir2* ...>
        only run tests from these sub-dirs or their children
        default = none
      -customonly <file1 file2* ...>
        only run tests from sub-dirs that contain these files
        default = none
      -custom <file_prefix>
        read options from this file_prefix plus test name in each sub-dir, if it exists
        valid options are: launch, descriptors, tolerance, error_norm, relative_error
        the number of launches and descriptors must match
        lines in file have syntax like:
          descriptors = ["1","4","8"]
          error_norm = "L1"
        default = "options"
      -error_norm <"L1" or "L2" or "max">
        metric for comparing a column of output to gold standard answer
        these are vector norms, treating the column as a vector
        default = "max"
      -relative_error <"True" or "False">
        treat tolerance as a relative error or not
        default = "False"
      -tolerance <float>
        column difference > tolerance = fail, else success
        default = 1.0e-7
      -logread <dir module>
        path for where to find the log-file reading Python module
        default = . log (log.py in this dir or somewhere Python can find it)
