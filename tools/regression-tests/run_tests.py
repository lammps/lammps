#!/usr/bin/env python3
'''
UPDATE: September 4, 2024:
  Launching the LAMMPS binary under testing using a configuration defined in a yaml file (e.g. config.yaml).
  Comparing the output thermo with that in the existing log file (with the same nprocs)
    + data in the log files are extracted and converted into yaml data structure
    + using the in place input scripts, no need to add REG markers to the input scripts

With the current features, users can:
    + specify which LAMMPS binary version to test (e.g., the version from a commit, or those from `lammps-testing`)
    + specify the examples subfolders (thus the reference log files) seperately (e.g. from other LAMMPS versions or commits)
    + specify the list of examples input scripts to test
    + specify tolerances for individual quantities for any input script to override the global values
    + launch tests with `mpirun` with all supported command line features (multiple procs, multiple paritions, and suffixes)
    + skip certain input files (whose names match specified patterns) if not interested, or packaged not installed, or no reference log file exists
    + set a timeout for every input script run if they may take too long
    + skip numerical checks if the goal is just to check if the runs do not fail

Some benefits include:

    + separating regression testing from building LAMMPS
    + performing quick and full regression tests
    + keeping track of the testing progress to resume the testing from the last checkpoint (skipping completed runs)
    + distributing the input list across multiple processes by
      splitting the list of input scripts into separate runs (there are ~800 input scripts under the top-level examples)
    + generating new reference log files if desirable

Input arguments:
    + the path to a LAMMPS binary (can be relative to the working directory)
    + a test configuration file (see tools/regression-tests/config.yaml for an example)
    + a text file that lists of folders where the input scripts reside and how many of them line by line, or
      a text file that list of input scripts, or
      the path to the top-level examples

Output:
    + failure.yaml : list of the failed runs and reasons
    + progress.yaml: full testing results of the tested input scripts with the status (completed, failed or skipped)
                     with error messages (for failed runs), and walltime (in seconds)
    + output.xml   :    testing results in the JUnit XML format
    + run.log      :       screen output and error of individual runs

Limitations:
    - input scripts use thermo style multi (e.g., examples/peptide) do not work with the expected thermo output format
    - input scripts that require partition runs (e.g. examples/neb) need a separate config file, e.g. args: "--partition 3x1"
    - testing accelerator packages (GPU, INTEL, KOKKOS, OPENMP) need separate config files, "args: -sf omp -pk omp 4"

The following Python packages need to be installed into an activated environment:

    python3 -m venv testing-env
    source testing-env/bin/activate
    pip install numpy pyyaml junit_xml

Example usage (aka, tests for this script):

    1) Simple use (using the provided tools/regression-tests/config.yaml and the examples/ folder at the top level)
           python3 run_tests.py --lmp-bin=build/lmp --config-file=tools/regression-tests/config.yaml

    2) Use a custom testing configuration
           python3 run_tests.py --lmp-bin=/path/to/lmp_binary --config-file=/path/to/config/file/config.yaml

    3) Specify a list of example folders
           python3 run_tests.py --lmp-bin=/path/to/lmp_binary --config-file=/path/to/config/file/config.yaml \
                --example-folders="/path/to/examples/melt;/path/to/examples/rigid"

       The example subfolders can also be loaded from a text file list_subfolders1.txt:
           python3 run_tests.py --lmp-bin=/path/to/lmp_binary --config-file=/path/to/config/file/config.yaml \
                --list-subfolders=list_subfolders1.txt --output-file=output1.txt --progress-file=progress1.yaml \
                --log-file=run1.log

    4) Specify a list of example input scripts (e.g. obtained from running tools/regression-tests/get-quick-list.py)
            python3 run_tests.py --lmp-bin=/path/to/lmp_binary --config-file=/path/to/config/file/config.yaml \
                --list-input=input_list.txt

    5) Test a LAMMPS binary with the whole top-level /examples folder in a LAMMPS source tree
           python3 run_tests.py --lmp-bin=/path/to/lmp_binary --examples-top-level=/path/to/lammps/examples
                --config-file=tools/regression-tests/config_serial.yaml

    6) Analyze the LAMMPS binary and whole top-level /examples folder in a LAMMPS source tree
       and generate separate input lists for 8 workers:
           python3 run_tests.py --lmp-bin=/path/to/lmp_binary --examples-top-level=/path/to/lammps/examples \
                --analyze --num-workers=8

       The output of this run is 8 files folder-list-[0-7].txt that lists the subfolders
       and 8 files input-list-[0-7].txt that lists the input scripts under the top-level example folders.
       With these lists, one can launch multiple instances of run_tests.py simultaneously
       each with a list of example subfolders (Case 3), or with a list of input scripts (Case 4).
'''

from argparse import ArgumentParser
import datetime
import fnmatch
import logging
import os
import random
import re
import subprocess
import sys
#from multiprocessing import Pool

# need "pip install numpy pyyaml"
import numpy as np
import yaml

# need "pip install junit_xml"
from junit_xml import TestSuite, TestCase

try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader

# infer top level LAMMPS dir from filename
LAMMPS_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', '..'))

# import git interface module
sys.path.append(os.path.realpath(os.path.join(LAMMPS_DIR, 'tools', 'regression-tests')))
import get_quick_list

'''
   data structure to store the test result
'''
class TestResult:
  def __init__(self, name, output=None, time=None, checks=0, status=None):
    self.name = name
    self.output = output
    self.time = time
    self.checks = 0
    self.status = status

'''
    Iterate over a list of input folders and scripts using the given lmp_binary and the testing configuration

    lmp_binary   : full path to the LAMMPS binary
    input_folder : the absolute path to the input files
    input_list   : list of the input scripts under the input_folder
    config       : the dict that contains the test configuration
    walltime_ref : reference walltime

    return
       results   : a list of TestResult objects
       stat      : a dictionary that lists the number of passed, skipped, failed tests
       progress_file: yaml file that stores the tested input script and status
       failure_file : file that reports the failed runs (a subset of progress_file)
       last_progress: the dictionary that shows the status of the last tests
       output_buf: placeholder for storing the output of a given worker
'''
def iterate(lmp_binary, input_folder, input_list, config, results, progress_file, failure_file, walltime_ref=1, verbose=False, last_progress=None, output_buf=None):

    num_tests = len(input_list)
    num_completed = 0
    num_passed = 0
    num_skipped = 0
    num_error = 0
    num_memleak = 0
    test_id = 0

    # using REG-commented input scripts, now turned off (False)
    using_markers = False
    EPSILON = np.float64(config['epsilon'])
    nugget = float(config['nugget'])
    genref = config['genref']
    compiler = config['compiler']
    use_valgrind = False
    if 'valgrind' in config['mpiexec']:
        use_valgrind = True

    # record all the failed runs
    failure = open(failure_file, "a")

    # iterate over the input scripts
    for input in input_list:

        # check if the progress file exists to append or create a new one
        if os.path.isfile(progress_file) == True:
            progress = open(progress_file, "a")
        else:
            progress = open(progress_file, "w")

        # walltime =   -2: skipped tests
        #              -1: failed tests
        #            >= 0: walltime in seconds (e.g. in.melt walltime = 0.2 seconds)
        walltime = -2

        # skip the input file if listed in the config file or matched with a pattern
        if 'skip' in config:
            if input in config['skip']:
                msg = "   + " + input + f" ({test_id+1}/{num_tests}): skipped as specified in {configFileName}"
                print(msg)
                logger.info(msg)
                progress.write(f"{input}: {{ folder: {input_folder}, status: \"skipped\", walltime: {walltime} }}\n")
                progress.close()
                num_skipped = num_skipped + 1
                test_id = test_id + 1
                continue

            matched_pattern = False
            for skipped_files in config['skip']:
                if '*' in skipped_files:
                    # check input script name e.g. in.*_imd*
                    if fnmatch.fnmatch(input, skipped_files):
                        matched_pattern = True
                        break

            if matched_pattern == True:
                msg = "   + " + input + f" ({test_id+1}/{num_tests}): skipped as specified in {configFileName}"
                print(msg)
                logger.info(msg)
                progress.write(f"{input}: {{ folder: {input_folder}, status: \"skipped\", walltime: {walltime} }}\n")
                progress.close()
                num_skipped = num_skipped + 1
                test_id = test_id + 1
                continue

        # also skip if the test already completed as marked in the progress file
        if input in last_progress:
            status = last_progress[input]['status']
            if 'completed' in status or 'numerical checks skipped' in status:
                msg = "  + " + input + f" ({test_id+1}/{num_tests}): marked as completed or numerical checks skipped (see {progress_file})"
                logger.info(msg)
                print(msg)
                # No need to write to progress again that the run is completed
                progress.close()
                num_skipped = num_skipped + 1
                test_id = test_id + 1
                continue

            if 'packaged not installed' in status:
                msg = "  + " + input + f" ({test_id+1}/{num_tests}): due to package not installed (see {progress_file})"
                logger.info(msg)
                print(msg)
                # No need to write to progress again that the run gets error due to missing packages
                progress.close()
                num_skipped = num_skipped + 1
                test_id = test_id + 1
                continue

        # if annotating input scripts with REG markers is True
        if using_markers == True:
            input_test = 'test.' + input
            if os.path.isfile(input) == True:
                if has_markers(input):
                    process_markers(input, input_test)

                else:
                    print(f"WARNING: {input} does not have REG markers")
                    input_markers = input + '.markers'
                    # if the input file with the REG markers does not exist
                    #   attempt to plug in the REG markers before each run command
                    if os.path.isfile(input_markers) == False:
                        cmd_str = "cp " + input + " " + input_markers
                        os.system(cmd_str)
                        generate_markers(input, input_markers)
                        process_markers(input_markers, input_test)

        else:
            # else the same file name for testing
            input_test = input

        str_t = "   + " + input_test + f" ({test_id+1}/{num_tests})"
        logger.info(str_t)
        print(str_t)

        # check if a reference log file exists in the current folder: log.DDMMMYY.basename.g++.[nprocs]
        # assuming that input file names start with "in." (except in.disp, in.disp2 and in.dos in phonon/)
        basename = input_test[3:]
        ref_logfile_exist = False

        # if there are multiple log files for different number of procs, pick the maximum number
        cmd_str = "ls log.*"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        logfile_list = p.stdout.split('\n')
        logfile_list.remove('')

        max_np = 1
        for file in logfile_list:
            # looks for pattern log.{date}.{basename}.{compiler}.{nprocs}: log.[date].min.box.[compiler]].* vs log.[date].min.[compiler].*
            # get the date from the log files
            date = file.split('.',2)[1]
            compiler = file.rsplit('.',2)[1]
            pattern = f'log.{date}.{basename}.{compiler}.*'
            if fnmatch.fnmatch(file, pattern):
                p = file.rsplit('.', 1)
                if p[1].isnumeric():
                    # if using valgrind or running in serial, then use the log file with 1 proc
                    if use_valgrind == True or config['mpiexec'] == "":
                        if int(p[1]) == 1:
                            max_np = int(p[1])
                            ref_logfile_exist = True
                            thermo_ref_file = file
                            break
                    else:
                        if max_np <= int(p[1]):
                            max_np = int(p[1])
                            ref_logfile_exist = True
                            thermo_ref_file = file

        # if there is no ref log file and not running with valgrind
        if ref_logfile_exist == False and use_valgrind == False:
            max_np = 4

        saved_nprocs = config['nprocs']

        # if the maximum number of procs is different from the value in the configuration file
        #      then override the setting for this particular input script
        if max_np != int(config['nprocs']):
            config['nprocs'] = str(max_np)

        # store the value of nprocs
        nprocs = int(config['nprocs'])

        # if valgrind is used for mem check, the run command will be
        #   mpirun -np 1 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes /path/to/lmp_binary -in in.txt
        # so both mpiexec_numproc_flag and nprocs are empty
        if use_valgrind == True:
            config['nprocs'] = ""
            config['mpiexec_numproc_flag'] = ""
            nprocs = 1

        # default walltime value of failed tests
        walltime = -1

        result = TestResult(name=input, output="", time="", status="passed")

        # run the LAMMPS binary with the input script
        cmd_str, output, error, returncode, logfilename = execute(lmp_binary, config, input_test)

        # restore the nprocs value in the configuration
        config['nprocs'] = saved_nprocs

        # check if the output contains ERROR
        #   there might not be a log file generated at this point, or only the log file contains only the date line
        if "ERROR" in output:
            error_line = ""
            for line in output.split('\n'):
                if "ERROR" in line:
                    error_line = line
                    break
            logger.info(f"     The run terminated with {input_test} gives the following output:")
            logger.info(f"       {error_line}")
            if "Unrecognized" in output:
                result.status = f"failed, unrecognized command, package not installed, {error_line}"
            elif "Unknown" in output:
                result.status = f"failed, unknown command, package not installed, {error_line}"
            else:
                result.status = f"failed, {error_line}."

            logger.info(f"     Output:")
            logger.info(f"     {output}")
            logger.info(f"     Failed with {input_test}.\n")
            num_error = num_error + 1

            results.append(result)
            print(f"{result.status}")

            msg = f"{input}: {{ folder: {input_folder}, status: \"{result.status}\", walltime: {walltime} }}\n"
            progress.write(msg)
            progress.close()
            failure.write(msg)

            test_id = test_id + 1
            continue

        # check if a log file log.{basename}.{nprocs} exists in the current folder
        if os.path.isfile(logfilename) == False:
            msg = f"    failed, no log.{basename}.{nprocs} generated with {input_test} with return code {returncode}.\n"
            print(msg)
            logger.info(msg)
            logger.info(f"    Output:")
            logger.info(f"    {output}")
            logger.info(f"    Error:\n{error}")

            msg = f"{input}: {{ folder: {input_folder}, status: \"failed, no log file generated\", walltime: {walltime} }}\n"
            progress.write(msg)
            progress.close()
            failure.write(msg)

            num_error = num_error + 1
            test_id = test_id + 1
            continue
        else:
            # generate a new log file whose name has the format of log.{date}.{basename}.{compiler}.{nprocs}
            if genref == True:
                dmy = datetime.datetime.now()
                date = dmy.strftime("%d%b%y")
                # assume g++ for now, but is be available from running "lmp_binary -h"
                compiler = "g++"
                cmd_str = f"cp log.{basename}.{nprocs} log.{date}.{basename}.{compiler}.{nprocs}"
                p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)

        # if skip numerical checks, then skip the rest
        if skip_numerical_check == True:
            msg = "completed, skipping numerical checks"
            if use_valgrind == True:
                if "All heap blocks were freed" in error:
                    msg += ", no memory leak"
                else:
                    msg += ", memory leaks detected"
                    num_memleak = num_memleak + 1
            result.status = msg
            results.append(result)

            msg = f"{input}: {{ folder: {input_folder}, status: \"{msg}\", walltime: {walltime} }}\n"
            progress.write(msg)
            progress.close()
            failure.write(msg)

            # count the number of completed runs
            num_completed = num_completed + 1
            test_id = test_id + 1
            continue

        # if there is no ERROR in the output, but there is no Total wall time printed out
        if "Total wall time" not in output:
            msg = f"     failed, no Total wall time in the output.\n"
            print(msg)
            logger.info(msg)
            logger.info(f"\n{input_test}:")
            logger.info(f"\n    Output:\n{output}")
            logger.info(f"\n    Error:\n{error}")

            msg = f"{input}: {{ folder: {input_folder}, status: \"failed, no Total wall time in the output, {error}\", walltime: {walltime} }}\n"
            progress.write(msg)
            progress.close()
            failure.write(msg)

            num_error = num_error + 1
            test_id = test_id + 1
            continue

        # NOTE: Total wall time could be 00:00:00 whereas Loop time is non-zero seconds
        walltime_norm = 1.0
        for line in output.split('\n'):
            if "Total wall time" in line:
                walltime_str = line.split('time:')[1]
                hms = walltime_str.split(':')
                hours = float(hms[0])
                minutes = float(hms[1])
                seconds = float(hms[2])
                walltime = hours * 3600.0 + minutes * 60.0 + seconds
                walltime_norm = float(walltime) / float(walltime_ref)
                break

        # if there is no Step or no Loop printed out
        if "Step" not in output or "Loop" not in output:
            msg = f"    completed, but no Step nor Loop in the output.\n"
            print(msg)
            logger.info(msg)
            logger.info(f"\n{input_test}:")
            logger.info(f"\n    Output:\n{output}")
            logger.info(f"\n    Error:\n{error}")

            msg = f"{input}: {{ folder: {input_folder}, status: \"completed, but no Step nor Loop in the output.\", walltime: {walltime}, walltime_norm: {walltime_norm} }}\n"
            progress.write(msg)
            progress.close()
            failure.write(msg)

            num_error = num_error + 1
            test_id = test_id + 1
            continue

        # parse thermo output in the log file from the run
        thermo = extract_data_to_yaml(logfilename)
        num_runs = len(thermo)

        # the run completed normally but the log file may not be friendly for parsing into YAML format
        if num_runs == 0:
            logger.info(f"     The run terminated with {input_test} gives the following output:")
            logger.info(f"     {output}")

            msg = "completed"
            if 'valgrind' in config['mpiexec']:
                if "All heap blocks were free" in error:
                    msg += ", no memory leak"
                else:
                    msg += ", memory leaks detected"
                    num_memleak = num_memleak + 1

            result.status = msg + f", error parsing {logfilename} into YAML"
            results.append(result)
            progress.write(f"{input}: {{ folder: {input_folder}, status: \"{result.status}\", walltime: {walltime}, walltime_norm: {walltime_norm} }}\n")
            progress.close()

            if verbose == True:
                print(result.status)

            num_completed = num_completed + 1
            test_id = test_id + 1
            continue

        # At this point, the run completed without trivial errors, proceed with numerical checks for thermo output
        # check if there is a reference log file for this input
        if ref_logfile_exist:
            # parse the thermo output in reference log file
            thermo_ref = extract_data_to_yaml(thermo_ref_file)
            if thermo_ref:
                num_runs_ref = len(thermo_ref)
            else:
                # thhe thermo_ref dictionary is empty
                logger.info(f"    failed, error parsing the reference log file {thermo_ref_file}.")
                result.status = "skipped numerical checks due to parsing the reference log file"
                results.append(result)
                progress.write(f"{input}: {{ folder: {input_folder}, status: \"completed, numerical checks skipped, unsupported log file format\", walltime: {walltime}, walltime_norm: {walltime_norm} }}\n")
                progress.close()
                num_completed = num_completed + 1
                num_error = num_error + 1
                test_id = test_id + 1
                continue
        else:
            msg = f"       failed, cannot find the reference log file for {input_test} with the expected format log.[date].{basename}.*.[nprocs]"
            logger.info(msg)
            print(msg)
            # attempt to read in the thermo yaml output from the working directory (the following section will be deprecated)
            thermo_ref_file = 'thermo.' + input + '.yaml'
            file_exist = os.path.isfile(thermo_ref_file)
            if file_exist == True:
                thermo_ref = extract_thermo(thermo_ref_file)
                num_runs_ref = len(thermo_ref)
            else:
                # most likely to reach here if the reference log file does not exist
                logger.info(f"       {thermo_ref_file} also does not exist in the working directory.")
                result.status = "skipped due to missing the reference log file"
                results.append(result)

                msg = f"{input}: {{ folder: {input_folder}, status: \"completed, numerical checks skipped due to missing the reference log file\", walltime: {walltime}, walltime_norm: {walltime_norm} }}\n"
                progress.write(msg)
                progress.close()
                failure.write(msg)
                num_completed = num_completed + 1
                num_error = num_error + 1
                test_id = test_id + 1
                continue

        logger.info(f"     Comparing thermo output from {logfilename} against the reference log file {thermo_ref_file}")

        # check if the number of runs matches with that in the reference log file
        # maybe due to some changes to the input where the ref log file is not updated yet
        if num_runs != num_runs_ref:
            logger.info(f"     ERROR: Number of runs in {logfilename} ({num_runs}) is different from that in the reference log ({num_runs_ref})."
                        " Check README in the folder, possibly due to using mpirun with partitions or parsing the wrong reference log file.")
            result.status = "failed, incomplete runs"
            results.append(result)
            progress.write(f"{input}: {{ folder: {input_folder}, status: \"{result.status}\", walltime: {walltime}, walltime_norm: {walltime_norm} }}\n")
            progress.close()
            num_error = num_error + 1
            test_id = test_id + 1
            continue

        # check if the number of fields match with that in the reference log file in the first run
        # due to some changes to the input where the ref log file is not updated yet
        # for early exit
        num_fields = len(thermo[0]['keywords'])
        num_fields_ref = len(thermo_ref[0]['keywords'])
        if num_fields != num_fields_ref:
            logger.info(f"     failed, number of thermo colums in {logfilename} ({num_fields}) is different from that in the reference log ({num_fields_ref}) in the first run.")
            logger.info(f"     Check both log files for more details.")
            result.status = "failed, mismatched columns in the log files"
            results.append(result)
            progress.write(f"{input}: {{ folder: {input_folder}, status: \"{result.status}\", walltime: {walltime}, walltime_norm: {walltime_norm} }}\n")
            progress.close()
            num_error = num_error + 1
            test_id = test_id + 1
            continue

        # comparing output vs reference values
        width = 20
        if verbose == True:
            print("        Quantities".ljust(width) + "Output".center(width) + "Reference".center(width) +
                "Abs Diff Check".center(width) +  "Rel Diff Check".center(width))

        # check if overrides for this input scipt is specified
        overrides = {}
        if 'overrides' in config:
            if input_test in config['overrides']:
                overrides = config['overrides'][input_test]

        # iterate through num_runs

        num_abs_failed = 0
        num_rel_failed = 0
        failed_abs_output = []
        failed_rel_output = []
        num_checks = 0
        mismatched_columns = False

        for irun in range(num_runs):
            num_fields = len(thermo[irun]['keywords'])
            num_fields_ref = len(thermo_ref[irun]['keywords'])
            if num_fields != num_fields_ref:
                logger.info(f"     failed: Number of thermo columns in {logfilename} ({num_fields})")
                logger.info(f"     is different from that in the reference log ({num_fields_ref}) in run {irun}.")
                mismatched_columns = True
                continue

            # get the total number of the thermo output lines
            nthermo_steps = len(thermo[irun]['data'])

            # get the output at the last timestep
            thermo_step = nthermo_steps - 1

            # iterate over the fields
            for i in range(num_fields):
                quantity = thermo[irun]['keywords'][i]

                val = thermo[irun]['data'][thermo_step][i]
                ref = thermo_ref[irun]['data'][thermo_step][i]
                abs_diff = abs(float(val) - float(ref))

                if abs(float(ref)) > EPSILON:
                    rel_diff = abs(float(val) - float(ref))/abs(float(ref))
                else:
                    rel_diff = abs(float(val) - float(ref))/abs(float(ref)+nugget)

                abs_diff_check = "PASSED"
                rel_diff_check = "PASSED"

                if quantity in config['tolerance'] or quantity in overrides:

                    if quantity in config['tolerance']:
                        abs_tol = float(config['tolerance'][quantity]['abs'])
                        rel_tol = float(config['tolerance'][quantity]['rel'])

                    # overrides the global tolerance values if specified
                    if quantity in overrides:
                        abs_tol = float(overrides[quantity]['abs'])
                        rel_tol = float(overrides[quantity]['rel'])

                    num_checks = num_checks + 2
                    if abs_diff > abs_tol:
                        abs_diff_check = "FAILED"
                        reason = f"Run {irun}: {quantity}: actual ({abs_diff:0.2e}) > expected ({abs_tol:0.2e})"
                        failed_abs_output.append(f"{reason}")
                        num_abs_failed = num_abs_failed + 1
                    if rel_diff > rel_tol:
                        rel_diff_check = "FAILED"
                        reason = f"Run {irun}: {quantity}: actual ({rel_diff:0.2e}) > expected ({rel_tol:0.2e})"
                        failed_rel_output.append(f"{reason}")
                        num_rel_failed = num_rel_failed + 1
                else:
                    # N/A means that tolerances are not defined in the config file
                    abs_diff_check = "N/A"
                    rel_diff_check = "N/A"

                if verbose == True and abs_diff_check != "N/A"  and rel_diff_check != "N/A":
                    print(f"        {thermo[irun]['keywords'][i].ljust(width)} {str(val).rjust(20)} {str(ref).rjust(20)} {abs_diff_check.rjust(20)} {rel_diff_check.rjust(20)}")

        # after all runs completed, or are interrupted in one of the runs (mismatched_columns = True)

        if mismatched_columns == True:
            msg = f"     mismatched log files after the first run. Check both log files for more details."
            print(msg)
            logger.info(msg)
            result.status = "thermo checks failed due to mismatched log files after the first run"

        result.status = ""
        if num_abs_failed > 0:
            msg = f"     {num_abs_failed} abs diff checks failed."
            print(msg)
            logger.info(msg)
            #result.status = f"abs_diff_failed: {num_abs_failed}, "
            if verbose == True:
                for out in failed_abs_output:
                    print(f"        - {out}")

        if num_rel_failed > 0:
            msg = f"     {num_rel_failed} rel diff checks failed."
            print(msg)
            logger.info(msg)
            #result.status += f"rel_diff_failed: {num_rel_failed}"
            if verbose == True:
                for out in failed_rel_output:
                    print(f"        - {out}")

        if num_abs_failed == 0 and num_rel_failed == 0:
            msg = f"     all {num_checks} checks passed."
            print(msg)
            logger.info(msg)
            #result.status = f"all {num_checks} checks passed."
            num_passed = num_passed + 1
        else:
            num_error = num_error + 1

        result.status = f"abs_diff_failed: {num_abs_failed}, rel_diff_failed: {num_rel_failed}"
        results.append(result)

        # check if memleak detects from valgrind run (need to replace "mpirun" -> valgrind --leak-check=yes mpirun")
        msg = "completed"
        if use_valgrind == True:
            if "All heap blocks were freed" in error:
                msg += ", no memory leak"
            else:
                msg += ", memory leaks detected"
                num_memleak = num_memleak + 1

        progress.write(f"{input}: {{ folder: {input_folder}, status: \"{msg}\", failed_checks: {{ {result.status} }}, walltime: {walltime}, walltime_norm: {walltime_norm} }}\n")
        progress.close()

        # write to failure if there is any numerical failed check
        if num_abs_failed > 0 or num_rel_failed > 0:
            failure.write(f"{input}: {{ folder: {input_folder}, status: \"{msg}\", failed_checks: {{ {result.status} }}, walltime: {walltime}, walltime_norm: {walltime_norm} }}\n")

        # count the number of completed runs
        num_completed = num_completed + 1
        test_id = test_id + 1

    # close the failure file
    failure.close()

    stat = { 'num_completed': num_completed,
             'num_passed': num_passed,
             'num_skipped': num_skipped,
             'num_error': num_error,
             'num_memleak':  num_memleak,
           }
    return stat

# HELPER FUNCTIONS
'''
  get the thermo output from a log file with thermo style yaml

  yamlFileName: input YAML file with thermo structured
      as described in https://docs.lammps.org/Howto_structured_data.html
  return: thermo, which is a list containing a dictionary for each run
      where the tag "keywords" maps to the list of thermo header strings
      and the tag data has a list of lists where the outer list represents the lines
      of output and the inner list the values of the columns matching the header keywords for that step.
'''
def extract_thermo(yamlFileName):
    docs = ""
    with open(yamlFileName) as f:
        for line in f:
            m = re.search(r"^(keywords:.*$|data:$|---$|\.\.\.$|  - \[.*\]$)", line)
            if m: docs += m.group(0) + '\n'
        thermo = list(yaml.load_all(docs, Loader=Loader))
        return thermo


'''
    Convert an existing log file into a thermo yaml style log
    inputFileName = a provided log file in an examples folder (e.g. examples/melt/log.8Apr21.melt.g++.4)
    return a YAML data structure as if loaded from a thermo yaml file
'''
def extract_data_to_yaml(inputFileName):
    with open(inputFileName, 'r') as file:
        data = file.read()
        lines = data.splitlines()
        reading = False
        data = []
        docs = ""
        num_thermo_cols = 0
        for line in lines:
            if "Step" in line and line[0] != '#':
                line.strip()
                keywords = line.split()
                num_thermo_cols = len(keywords)
                reading = True
                docs += "---\n"
                docs += str("keywords: [")
                for word in enumerate(keywords):
                  docs += "'" + word[1] + "', "
                docs += "]\n"
                docs += "data:\n"
            if "Loop" in line:
                reading = False
                docs += "...\n"

            if reading == True and "Step" not in line:
                if  "WARNING" in line:
                    continue
                data = line.split()
                if len(data) != num_thermo_cols:
                    continue
                docs += " - ["
                for field in enumerate(data):
                    docs += field[1] + ", "
                docs += "]\n"

    # load the docs into a YAML data struture
    #print(docs)
    thermo = {}
    try:
        yaml_struct = yaml.load_all(docs, Loader=Loader)
        thermo = list(yaml_struct)
    except yaml.YAMLError as exc:
        if hasattr(exc, 'problem_mark'):
            mark = exc.problem_mark
            msg = f"    Error parsing {inputFileName} at line {mark.line}, column {mark.column+1}."
            print(msg)
            logger.info(msg)
            logger.info(docs)
            return thermo
        else:
            msg = f"    Something went wrong while parsing {inputFileName}."
            print(msg)
            logger.info(msg)
            logger.info(docs)
            return thermo
    return thermo

'''
    return a dictionary of the list of installed packages, OS, GitInfo, compiler and compile_flags
'''
def get_lammps_build_configuration(lmp_binary):
    cmd_str = lmp_binary + " -h"
    p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
    output = p.stdout.split('\n')
    packages = ""
    reading = False
    operating_system = ""
    GitInfo = ""
    compiler = "g++"
    compiler_full = ""
    row = 0
    for line in output:
        if line != "":
            if line == "Installed packages:":
                reading = True
                n = row
            if "List of individual style options" in line:
                reading = False
            if reading == True and row > n:
                packages += line.strip() + " "

        if "OS:" in line:
            operating_system = line
        if "Git info" in line:
            GitInfo = line
        if "Compiler" in line:
            compiler_full = line
            if "GNU" in line:
                compiler = "g++"
            if "Intel" in line: 
                compiler = "icc"
        row += 1

    packages = packages.strip()

    row = 0
    compile_flags = ""
    for line in output:
        if line != "":
            if "-DLAMMPS" in line:
                compile_flags += " " + line.strip()

        row += 1

    installed_packages = packages.split(" ")
    build_config = {
        'installed_packages': installed_packages,
        'operating_system': operating_system,
        'git_info': GitInfo, 
        'compiler': compiler,
        'compiler_full': compiler_full,
        'compile_flags': compile_flags,
    }

    return build_config

'''
    launch LAMMPS using the configuration defined in the dictionary config with an input file
    return
       - cmd_str:     the complete command used to launch LAMMPS with the input
       - stdout:      stdout of the process
       - stderr:      stderr of the process
       - errorcode:   error code returned by the process
       - logfilename: the log file name for the given input
                      to avoid duplicate writes to log.lammps if multiple workers execute in the same folder
'''
def execute(lmp_binary, config, input_file_name, generate_ref=False):
    cmd_str = ""
    # check if mpiexec/mpirun is used
    if config['mpiexec']:
        cmd_str += config['mpiexec'] + " " + config['mpiexec_numproc_flag'] + " " + config['nprocs'] + " "

    # write to a log file with format log.{basename}.{nprocs}
    basename = input_file_name[3:]
    logfilename = f"log.{basename}.{config['nprocs']}"

    cmd_str += lmp_binary + " -in " + input_file_name + " " + config['args'] + " -log " + logfilename

    logger.info(f"     Executing: {cmd_str}")
    # set a timeout (in seconds) for each run
    timeout = 60
    if 'timeout' in config:
        if config['timeout'] != "":
            timeout = int(config['timeout'])

    try:
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True, timeout=timeout)
        return cmd_str, p.stdout, p.stderr, p.returncode, logfilename

    except subprocess.TimeoutExpired:
        msg = f"     Timeout for: {cmd_str} ({timeout}s expired)"
        logger.info(msg)
        print(msg)

    error_str = f"timeout ({timeout}s expired)"
    return cmd_str, "", error_str, -1, logfilename

'''
   get the reference walltime by running the lmp_binary with config with an input script in the bench/ folder
      in.lj is suitable as it doesn't need any potential file, nor any extra packages
'''
def get_reference_walltime(lmp_binary, config):
    cmd_str = ""
    # check if mpiexec/mpirun is used
    if config['mpiexec']:
        cmd_str += config['mpiexec'] + " " + config['mpiexec_numproc_flag'] + " " + config['nprocs'] + " "

    # guess the build folder path
    lmp_build_folder = lmp_binary.rsplit('/', 1)[0]

    # guess the bench folder
    lmp_bench_folder = lmp_build_folder + "/../bench/"

    # run with replicate for a copple of seconds long run
    cmd_str += lmp_binary + " -in " + lmp_bench_folder + "in.lj -v x 2 -v y 2 -v z 1 " + config['args']

    logger.info(f"     Executing for reference walltime: {cmd_str}")

    # walltime = -1 indicates some timeout (issues)
    walltime = -1

    # set a timeout for this reference run
    timeout = 60
    output = ""
    try:
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True, timeout=timeout)
        output = p.stdout

    except subprocess.TimeoutExpired:
        msg = f"     Timeout for: {cmd_str} ({timeout}s expired)"
        logger.info(msg)
        print(msg)

    for line in output.split('\n'):
        if "Total wall time" in line:
            walltime_str = line.split('time:')[1]
            hms = walltime_str.split(':')
            hours = float(hms[0])
            minutes = float(hms[1])
            seconds = float(hms[2])
            walltime = hours * 3600.0 + minutes * 60.0 + seconds

    logger.info(f"     Reference walltime, sec = {walltime}")

    return walltime

'''
    infer the tools/regression-tests folder from the absolute path to lmp_binary
    return the default config file path tools/regression-tests/config.yaml
'''
def get_default_config(lmp_binary):
    # guess the build folder path
    lmp_build_folder = lmp_binary.rsplit('/', 1)[0]

    # guess the tools/regression-tests folder
    regression_tests_folder = lmp_build_folder + "/../tools/regression-tests/"

    defaultConfigFile = regression_tests_folder + "config.yaml"
    return defaultConfigFile

'''
    split a list into a list of N sublists

    NOTE:
    To map a function to individual workers with multiprocessing.Pool:

    def func(input1, input2, output_buf):
        # do smth
        return result

    # args is a list of num_workers tuples, each tuple contains the arguments passed to the function executed by a worker
    args = []
    for i in range(num_workers):
        args.append((input1, input2, output_buf))

    with Pool(num_workers) as pool:
        results = pool.starmap(func, args)
'''
def divide_into_N(original_list, N):
    size = np.ceil(len(original_list) / N)
    b = []
    for i in range(0, N):
        start = int(i * size)
        end = int(start + size)
        l = original_list[start:end]
        b.append(l)
    return b

'''
    process the #REG markers in an input script, add/replace with what follows each marker

    inputFileName:  LAMMPS input file with comments #REG:ADD and #REG:SUB as markers
    outputFileName: modified input file ready for testing
'''
def process_markers(inputFileName, outputFileName):
  # read in the script
  with open(inputFileName, 'r') as file:
    data = file.read()

    # replace #REG:ADD with empty string (i.e. adding the text at the end of the line)
    data = data.replace("#REG:ADD", "")

    # replace the line contaning #REG:SUB with a line with the text that follows this marker
    data = data.splitlines()
    separator="#REG:SUB"
    out = []
    for line in data:
      s = line.split(separator)
      if len(s) < 2:
        out.append(line)
      else:
        out.append(s[1])

  # write data to the new script
  with open(outputFileName, 'w') as file:
    for line in out:
      file.write(line + "\n")


'''
    attempt to insert the #REG markers before each run command
    #REG:ADD thermo 10
    #REG:ADD thermo_style yaml

    inputFileName:  provided LAMMPS input file
    outputFileName: modified input file ready for testing
'''
def generate_markers(inputFileName, outputFileName):
    # read in the script
    with open(inputFileName, 'r') as file:
        data = file.read()

        lines = data.splitlines()
        out = []
        for line in lines:
            s = line.split()
            if len(s) > 0:
                if s[0] == "run":
                    out.append("    #REG:ADD thermo 10")
                    out.append("    #REG:ADD thermo_style yaml")
            out.append(line)

    # write data to the new script
    with open(outputFileName, 'w') as file:
        for line in out:
          file.write(line + "\n")

'''
    check if any input script has any #REG markers
'''
def has_markers(inputFileName):
    with open(inputFileName) as f:
        if '#REG' in f.read():
          return True
    return False


'''
    Main entry
'''
if __name__ == "__main__":

    # default values
    lmp_binary = ""
    configFileName = "config.yaml"
    example_subfolders = []
    example_inputs = []
    example_toplevel = ""
    genref = False
    verbose = False
    output_file = "output.xml"
    progress_file = "progress.yaml"
    failure_file = "failure.yaml"
    log_file = "run.log"
    list_input = ""
    list_subfolders = ""
    analyze = False
    quick = False
    quick_branch = "origin/develop"
    quick_max = 50
    quick_reference = os.path.join(LAMMPS_DIR, 'tools', 'regression-tests', 'reference.yaml')

    # distribute the total number of input scripts over the workers
    num_workers = 1

    # parse the arguments
    parser = ArgumentParser()
    parser.add_argument("--lmp-bin", dest="lmp_binary", default="", help="LAMMPS binary")
    parser.add_argument("--config-file", dest="config_file", default="", help="Configuration YAML file")
    parser.add_argument("--examples-top-level", dest="example_toplevel", default="", help="Examples top-level")
    parser.add_argument("--example-folders", dest="example_folders", default="", help="Example subfolders")
    parser.add_argument("--list-input", dest="list_input", default="", help="File that lists the input scripts")
    parser.add_argument("--list-subfolders", dest="list_subfolders", default="", help="File that lists the subfolders")
    parser.add_argument("--num-workers", dest="num_workers", default=1, help="Number of workers")
    parser.add_argument("--output-file",dest="output", default=output_file, help="Output file")
    parser.add_argument("--log-file",dest="logfile", default=log_file, help="Log file")
    parser.add_argument("--progress-file",dest="progress_file", default=progress_file, help="Progress file")
    parser.add_argument("--failure-file",dest="failure_file", default=failure_file, help="Failure file")
    analyze = parser.add_mutually_exclusive_group()
    analyze.add_argument("--analyze",dest="analyze", action='store_true', default=False,
                        help="Analyze the testing folders and report statistics, not running the tests")
    analyze.add_argument("--quick", dest="quick", action='store_true', default=False,
                        help="Determine which test inputs have commands changed between a branch and the head")
    parser.add_argument("--quick-branch", dest="quick_branch", default=quick_branch,
                        help="Branch to which compare the current head to for changed styles")
    parser.add_argument("--quick-max", dest="quick_max", default=50,
                        help="Maximum number of inputs to randomly select")
    parser.add_argument("--quick-reference", dest="quick_reference", default=quick_reference,
                        help="Reference YAML file with progress data from full regression test run")
    parser.add_argument("--skip-numerical-check",dest="skip_numerical_check", action='store_true', default=False,
                        help="Skip numerical checks")
    parser.add_argument("--gen-ref",dest="genref", action='store_true', default=False,
                        help="Generating reference log files")
    parser.add_argument("--verbose",dest="verbose", action='store_true', default=False,
                        help="Verbose screen output")
    parser.add_argument("--resume",dest="resume", action='store_true', default=False,
                        help="Resume the test run from the list of inputs given the progress in progress.yaml")

    args = parser.parse_args()

    lmp_binary = os.path.abspath(args.lmp_binary)
    if len(args.config_file) > 0:
        configFileName = args.config_file
    else:
        configFileName = get_default_config(lmp_binary)

    output_file = args.output
    if int(args.num_workers) > 0:
        num_workers = int(args.num_workers)
    list_input = args.list_input
    list_subfolders = args.list_subfolders

    # example_toplevel is where all the examples subfolders reside
    if args.example_toplevel != "":
       example_toplevel = args.example_toplevel
    if args.example_folders != "":
        example_subfolders = args.example_folders.split(';')

    genref = args.genref
    verbose = args.verbose
    log_file = args.logfile
    analyze = args.analyze
    quick = args.quick
    quick_branch = args.quick_branch
    quick_max = int(args.quick_max)
    quick_reference = args.quick_reference
    skip_numerical_check = args.skip_numerical_check
    resume = args.resume
    progress_file = args.progress_file
    failure_file = args.failure_file

    # logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode="w")

    if len(example_subfolders) > 0:
        print("\nExample folders to test:")
        print(*example_subfolders, sep='\n')
    if example_toplevel != "":
        print("\nTop-level example folder:")
        print(f"  {example_toplevel}")
    if list_input != "":
        print("\nInput scripts to test as listed in the file:")
        print(f"  {list_input}")

    # Using in place input scripts
    inplace_input = True
    test_cases = []

    # generate list of input scripts with commands that have been changed
    if quick:
        headers = get_quick_list.changed_files_from_git(quick_branch)
        styles = get_quick_list.get_command_from_header(headers, LAMMPS_DIR)
        regex = get_quick_list.make_regex(styles)
        if regex:
            if not example_toplevel: example_toplevel = os.path.join(LAMMPS_DIR, 'examples')
            input_list = get_quick_list.get_examples_using_styles(regex, example_toplevel)
            msg = f"\nThere are {len(input_list)} input scripts with changed styles relative to branch {quick_branch}."
            msg += "\nChanged styles: " + str(styles)

            # read in refrence data from a previous test run
            with open(quick_reference, 'r') as f:
                reference = yaml.load(f, Loader=Loader)
            f.close()

            # trim previously failing run and runs that would take too long
            new_list = []
            keys = reference.keys()
            msg += "\nTrimming inputs using reference data from " + str(len(keys)) + " previous runs: "
            for infile in input_list:
                input = os.path.split(infile)[1]
                if input in keys:
                    if (reference[input]['walltime'] < 0.0):
                        # print("Skipping ", input, " for previous failure")
                        pass
                    elif (reference[input]['walltime'] > 29.0):
                        # print("Skipping ", input, " for wall time limit")
                        pass
                    else:
                        new_list.append(infile)
                else:
                    new_list.append(infile)
            input_list = new_list
            msg += "trimmed list has " + str(len(input_list)) + " entries"

            if len(input_list) > quick_max:
                input_list = random.sample(input_list, quick_max)
                msg += "\nTesting " + str(quick_max) + " randomly selected inputs"

            print(msg)
            logger.info(msg)

            # divide the list of input scripts into num_workers chunks
            sublists = divide_into_N(input_list, num_workers)

            # write each chunk to a file
            idx = 0
            for list_input in sublists:
                filename = f"input-list-{idx}.txt"
                with open(filename, "w") as f:
                    for inp in list_input:
                        f.write(inp + '\n')
                    f.close()
                idx = idx + 1
        else:
            msg = f"\nThere are no input scripts with changed styles relative to branch {quick_branch}."
            print(msg)
            logger.info(msg)
            for idx in range(0, num_workers):
                try:
                    os.remove(f"folder-list-{idx}.txt")
                except:
                    pass
                try:
                    os.remove(f"input-list-{idx}.txt")
                except:
                    pass
                filename = f"run-{idx}.log"
                with open(filename, "w") as f:
                    f.write('\n')
                f.close()
                filename = f"progress-{idx}.yaml"
                with open(filename, "w") as f:
                    f.write('\n')
                f.close()
                filename = f"output-{idx}.xml"
                with open(filename, "w") as f:
                    f.write('\n')
                f.close()
                filename = f"failure-{idx}.yaml"
                with open(filename, "w") as f:
                    f.write('\n')
                f.close()
        quit()

    # if the example folders are not specified from the command-line argument --example-folders
    # then use the path from --example-top-folder, or from the input-list read from a text file
    elif len(example_subfolders) == 0:

        # if the top level is specified
        if len(example_toplevel) != 0:
            # getting the list of all the input files because there are subfolders (e.g. PACKAGES) under the top level
            cmd_str = f"find {example_toplevel} -name \"in.*\" "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list.remove("")
            msg = f"\nThere are {len(input_list)} input scripts in total under the {example_toplevel} folder."
            print(msg)
            logger.info(msg)

            # get the input file list
            # TODO: generate a list of tuples, each tuple contains a folder list for a worker,
            #       then use multiprocessing.Pool starmap()
            folder_list = []
            for input in input_list:
                folder = input.rsplit('/', 1)[0]
                # unique folders in the list
                if folder not in folder_list:
                    folder_list.append(folder)

            # divide the list of folders into num_workers chunks
            sublists = divide_into_N(folder_list, num_workers)

            # write each chunk to a file
            idx = 0
            for list_input in sublists:
                filename = f"folder-list-{idx}.txt"
                with open(filename, "w") as f:
                    for folder in list_input:
                        # count the number of input scripts in each folder
                        cmd_str = f"ls {folder}/in.* | wc -l"
                        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
                        num_input = p.stdout.split('\n')[0]
                        f.write(folder + ' ' + num_input + '\n')
                    f.close()
                idx = idx + 1

            # working on all the folders for now
            example_subfolders = folder_list

            # divide the list of input scripts into num_workers chunks
            sublists = divide_into_N(input_list, num_workers)

            # write each chunk to a file
            idx = 0
            for list_input in sublists:
                filename = f"input-list-{idx}.txt"
                with open(filename, "w") as f:
                    for inp in list_input:
                        f.write(inp + '\n')
                    f.close()
                idx = idx + 1

        # if a list of subfolders is provided from a text file (list_subfolders from the command-line argument)
        elif len(list_subfolders) != 0:
            num_inputscripts = 0
            with open(list_subfolders, "r") as f:
                all_subfolders = f.read().splitlines()
                f.close()
                for line in all_subfolders:
                    if len(line) > 0:
                        # skip subfolders
                        if line[0] == '#':
                            continue
                        folder = line.split()[0]
                        example_subfolders.append(folder)
                        num_inputscripts += int(line.split()[1])
            msg = f"\nThere are {len(example_subfolders)} folders with {num_inputscripts} input scripts in total listed in {list_input}."
            print(msg)
            logger.info(msg)

        # if a list of input scripts is provided from a text file (list_input from the command-line argument)
        elif len(list_input) != 0:
            num_inputscripts = 0
            folder_list = []
            with open(list_input, "r") as f:
                all_inputs = f.read().splitlines()
                f.close()

                for line in all_inputs:
                    if len(line) > 0:
                        # skip input scripts
                        if line[0] == '#':
                            continue
                        input = line.split()[0]
                        folder = input.rsplit('/', 1)[0]
                        # unique folders in the list
                        if folder not in folder_list:
                            folder_list.append(folder)
                        example_inputs.append(input)
                        num_inputscripts += 1

            example_subfolders = folder_list
            msg = f"\nThere are {num_inputscripts} input scripts listed in {list_input}."
            print(msg)
            logger.info(msg)

        else:
            inplace_input = False

    # if analyze the example folders (and split into separate lists for top-level examples), not running any test
    if analyze == True:
        quit()

    # read in the configuration of the tests
    with open(configFileName, 'r') as f:
        config = yaml.load(f, Loader=Loader)
        absolute_path = os.path.abspath(configFileName)
        print(f"\nRegression test configuration file:\n  {absolute_path}")
        f.close()

    # check if lmp_binary is specified in the config yaml
    if lmp_binary == "":
        if config['lmp_binary'] == "":
            print("Needs a valid LAMMPS binary")
            quit()
        else:
            lmp_binary = os.path.abspath(config['lmp_binary'])

    # print out the binary info
    build_config = get_lammps_build_configuration(lmp_binary)
    packages = build_config['installed_packages']
    operating_system = build_config['operating_system']
    GitInfo = build_config['git_info']
    compiler = build_config['compiler']
    compiler_full = build_config['compiler_full']
    compile_flags = build_config['compile_flags']
    
    print("\nLAMMPS build info:")
    print(f"  - {operating_system}")
    print(f"  - {GitInfo}")
    print(f"  - {compiler_full}")
    print(f"  - Active compile flags: {compile_flags}")
    print(f"  - List of {len(packages)} installed packages:")
    all_pkgs = ""
    for p in packages:
        all_pkgs += p + " "
    print(all_pkgs)

    # augment config with additional keys
    config['compiler'] = compiler
    config['genref'] = genref

    all_results = []

    # save current working dir
    p = subprocess.run("pwd", shell=True, text=True, capture_output=True)
    pwd = p.stdout.split('\n')[0]
    pwd = os.path.abspath(pwd)
    print("\nWorking directory: " + pwd)

    progress_file_abs = pwd + "/" + progress_file
    last_progress = {}
    if resume == False:
        progress = open(progress_file_abs, "w")
        progress.close()
    else:
        try:
            progress = open(progress_file_abs, "r")
            last_progress = yaml.load(progress, Loader=Loader)
            progress.close()
        except Exception:
            print(f"    Cannot open progress file {progress_file_abs} to resume, rerun all the tests")

    # get a reference walltime
    walltime_ref = get_reference_walltime(lmp_binary, config)

    # record all the failure cases (overwrite if the file exists)
    failure_file_abs = pwd + "/" + failure_file
    failure = open(failure_file_abs, "w")
    failure.close()

    # initialize all the counters
    total_tests = 0
    completed_tests = 0
    passed_tests = 0
    skipped_tests = 0
    error_tests = 0
    memleak_tests = 0

    # default setting is to use inplace_input
    if inplace_input == True:

        # change dir to a folder under examples/
        # TODO: loop through the subfolders under examples/, depending on the installed packages

        '''
        args = []
        for i in range(num_workers):
            args.append((input1, input2, output))

        with Pool(num_workers) as pool:
            results = pool.starmap(func, args)
        '''

        for directory in example_subfolders:

            # change to the directory where the input script and data files are located
            print("-"*80)
            print("Entering " + directory)
            logger.info("Entering " + directory)
            os.chdir(directory)

            cmd_str = "ls in.*"
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            all_input_list = p.stdout.split('\n')
            all_input_list.remove('')

            # if the list of example input scripts is provided
            #   if an input script is not in the list, then remove it from input_list
            input_list = []
            if len(example_inputs) > 0:
                for inp in all_input_list:
                    full_path = directory + "/" + inp
                    if full_path in example_inputs:
                        input_list.append(inp)
            else:
                input_list = all_input_list

            print(f"{len(input_list)} input script(s) to be tested: {input_list}")
            total_tests += len(input_list)

            # iterate through the input scripts
            results = []
            stat = iterate(lmp_binary, directory, input_list, config,
                           results, progress_file_abs, failure_file_abs, walltime_ref, verbose, last_progress)

            completed_tests += stat['num_completed']
            skipped_tests += stat['num_skipped']
            passed_tests += stat['num_passed']
            error_tests += stat['num_error']
            memleak_tests += stat['num_memleak']

            # append the results to the all_results list
            all_results.extend(results)

            # get back to the working dir
            os.chdir(pwd)

    else:
        # or using the input scripts in the working directory -- for debugging purposes
        input_list=['in.lj']
        total_tests = len(input_list)
        results = []
        stat = iterate(lmp_binary, pwd, input_list, config, results, progress_file_abs)

        completed_tests = stat['num_completed']
        skipped_tests = stat['num_skipped']
        passed_tests = stat['num_passed']
        error_tests = stat['num_error']
        memleak_tests = stat['num_memleak']

        all_results.extend(results)

    # print out summary
    msg = "\nSummary:\n"
    msg += f"  Total number of input scripts: {total_tests}\n"
    msg += f"  - Skipped  : {skipped_tests}\n"
    msg += f"  - Failed   : {error_tests}\n"
    msg += f"  - Completed: {completed_tests}\n"

    # print notice to GitHub
    if 'GITHUB_STEP_SUMMARY' in os.environ:
        with open(os.environ.get('GITHUB_STEP_SUMMARY'), 'w') as f:
            print(f"Skipped: {skipped_tests}  Failed: {error_tests}  Completed: {completed_tests}", file=f)

    if memleak_tests < completed_tests and 'valgrind' in config['mpiexec']:
        msg += f"    - memory leak detected  : {memleak_tests}\n"
    if passed_tests <= completed_tests:
        msg += f"    - numerical tests passed: {passed_tests}\n"
    msg += "\nOutput:\n"
    msg += f"  - List of failed inputs         : {failure_file}\n"
    msg += f"  - Status of the tested inputs   : {progress_file}\n"
    msg += f"  - Running log with screen output: {log_file}\n"
    msg += f"  - Testing result in JUnit XML   : {output_file}\n"

    print(msg)

    # optional: need to check if junit_xml packaged is already installed in the env
    #   generate a JUnit XML file
    with open(output_file, 'w') as f:
        test_cases = []
        for result in all_results:
            #print(f"{result.name}: {result.status}")
            case = TestCase(name=result.name, classname=result.name)
            if result.status == "failed":
                case.add_failure_info(message="Actual values did not match expected ones.")
            if result.status == "skipped":
                case.add_skipped_info(message="Test was skipped.")
            if result.status == "error":
                case.add_skipped_info(message="Test run had errors.")
            test_cases.append(case)

        current_timestamp = datetime.datetime.now()
        ts = TestSuite(f"{configFileName}", test_cases, timestamp=current_timestamp)
        TestSuite.to_file(f, [ts], prettyprint=True)
