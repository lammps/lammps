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
    + skip certain input files (whose names match specified patterns) if not interested, or package not installed, or no reference log file exists
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
    + failure.yaml : list of the failed runs and reasons (a subset of progress.yaml)
    + progress.yaml: full testing results of the tested input scripts with the status (completed, failed or skipped)
                     with error messages (for failed runs), and walltime (in seconds)
    + output.xml   :    testing results in the JUnit XML format
    + run.log      :       screen output and error of individual runs

Both YAML files are written incrementally (one line per tested input) and parse as
a single YAML mapping keyed by the input script name. Multiple JUnit XML files from
distributed runs can be merged and summarized with merge_results.py in this folder.

Limitations:
    - input scripts use thermo style multi (e.g., examples/peptide) do not work with the expected thermo output format
    - input scripts that require partition runs (e.g. examples/neb) need a separate config file, e.g. args: "--partition 3x1"
    - testing accelerator packages (GPU, INTEL, KOKKOS, OPENMP) need separate config files, "args: -sf omp -pk omp 4"

The following Python packages need to be installed into an activated environment:

    python3 -m venv testing-env
    source testing-env/bin/activate
    pip install numpy pyyaml

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
import glob
import logging
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import xml.etree.ElementTree as ET
#from multiprocessing import Pool

# need "pip install numpy pyyaml"
import numpy as np
import yaml

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
   data structure to store the result of a single test run

   name    : name of the input script (e.g. in.melt)
   folder  : full path of the folder containing the input script
   category: normalized outcome for downstream reporting, one of
             'passed'    all numerical checks passed
             'failed'    the run completed but numerical checks failed
             'error'     the run did not complete (crash, timeout, missing log file)
             'completed' the run completed but no numerical checks were possible
                         (e.g. no reference log file, skipped numerical checks)
             'skipped'   the input script was not run at all
   status  : human readable description of the outcome (stored in the progress file)
   message : additional details (e.g. the individual failed checks)
   time    : walltime of the run in seconds (-1: failed, -2: skipped)
   checks  : number of performed numerical checks
   timeout : whether the run was killed after exceeding the timeout
   memleak : whether valgrind detected memory leaks
'''
class TestResult:
  def __init__(self, name, folder=None, category='skipped', status="", message="",
               time=-2.0, checks=0):
    self.name = name
    self.folder = folder
    self.category = category
    self.status = status
    self.message = message
    self.time = time
    self.checks = checks
    self.timeout = False
    self.memleak = False

'''
    condense a (possibly multi-line) text into a single line of limited length
    suitable for status messages in the progress file and JUnit XML attributes
'''
def shorten(text, maxlen=200):
    text = ' '.join(str(text).split())
    if len(text) > maxlen:
        text = text[:maxlen] + " ..."
    return text

'''
    write the results of the completed test runs into a JUnit XML file
    for downstream reporting (merging, summaries, websites)

    output_file: name of the JUnit XML file
    results    : list of TestResult objects
    suite_name : name of the test suite (e.g. the test configuration file)
    properties : dictionary with build and test configuration metadata
'''
def write_junit_xml(output_file, results, suite_name, properties=None):
    testsuite = ET.Element('testsuite')
    testsuite.set('name', suite_name)
    testsuite.set('timestamp', datetime.datetime.now().isoformat(timespec='seconds'))
    num_failures = 0
    num_errors = 0
    num_skipped = 0
    total_time = 0.0

    if properties:
        props = ET.SubElement(testsuite, 'properties')
        for key, value in properties.items():
            prop = ET.SubElement(props, 'property')
            prop.set('name', str(key))
            prop.set('value', shorten(value, maxlen=1000))

    for result in results:
        case = ET.SubElement(testsuite, 'testcase')
        case.set('name', result.name)
        # report the folder relative to the examples top-level folder, so that
        # results from different machines or checkouts can be compared
        folder = os.path.abspath(result.folder) if result.folder else ''
        marker = os.sep + 'examples' + os.sep
        if marker in folder:
            classname = folder.split(marker, 1)[1]
        else:
            classname = os.path.basename(folder)
        case.set('classname', classname.replace(os.sep, '/'))
        case.set('time', f"{max(float(result.time), 0.0):.3f}")
        total_time += max(float(result.time), 0.0)

        if result.category == 'failed':
            elem = ET.SubElement(case, 'failure')
            num_failures += 1
        elif result.category == 'error':
            elem = ET.SubElement(case, 'error')
            num_errors += 1
        elif result.category in ('skipped', 'completed'):
            # 'completed' means the run finished but nothing could be verified,
            # so it is reported as skipped rather than as passed
            elem = ET.SubElement(case, 'skipped')
            num_skipped += 1
        else:
            elem = None
        if elem is not None:
            elem.set('message', shorten(result.status))
            if result.message:
                elem.text = result.message

    testsuite.set('tests', str(len(results)))
    testsuite.set('failures', str(num_failures))
    testsuite.set('errors', str(num_errors))
    testsuite.set('skipped', str(num_skipped))
    testsuite.set('time', f"{total_time:.3f}")

    tree = ET.ElementTree(testsuite)
    ET.indent(tree)
    tree.write(output_file, encoding='UTF-8', xml_declaration=True)

'''
    Iterate over a list of input scripts using the given lmp_binary and the testing configuration

    lmp_binary   : full path to the LAMMPS binary
    input_folder : the absolute path to the input files
    input_list   : list of the input scripts under the input_folder
    config       : the dict that contains the test configuration
    results      : list of TestResult objects, one gets appended per input script
    progress_file: yaml file that stores the tested input script and status
    failure_file : file that reports the failed runs (a subset of progress_file)
    walltime_ref : reference walltime
    last_progress: the dictionary that shows the status of the last tests
    output_buf   : placeholder for storing the output of a given worker

    return
       stat      : a dictionary that lists the number of passed, skipped, failed tests
'''
def iterate(lmp_binary, input_folder, input_list, config, results, progress_file, failure_file, walltime_ref=1, verbose=False, last_progress=None, output_buf=None):

    num_tests = len(input_list)
    num_results_initial = len(results)
    test_id = 0

    EPSILON = np.float64(config['epsilon'])
    nugget = float(config['nugget'])
    genref = config['genref']
    compiler = config['compiler']
    use_valgrind = 'valgrind' in config['mpiexec']

    # record the outcome of a test: append it to the results list, write one line
    # to the progress file and, for failed or errored tests, to the failure file;
    # each line is a flow-style YAML mapping so that the whole file parses as a
    # single YAML mapping keyed by the input script names
    def record(result, walltime_norm=None, failed_checks=None, write_progress=True):
        results.append(result)
        if not write_progress:
            return
        entry = { 'folder': result.folder, 'status': result.status }
        if failed_checks is not None:
            entry['failed_checks'] = failed_checks
        entry['walltime'] = float(result.time)
        if walltime_norm is not None:
            entry['walltime_norm'] = float(walltime_norm)
        value = yaml.safe_dump(entry, default_flow_style=True,
                               sort_keys=False, width=1000000).strip()
        line = f"{result.name}: {value}\n"
        with open(progress_file, "a") as progress:
            progress.write(line)
        if result.category in ('failed', 'error'):
            with open(failure_file, "a") as failure:
                failure.write(line)

    # iterate over the input scripts
    for input in input_list:

        input_test = input
        result = TestResult(name=input, folder=input_folder)

        # skip the input file if listed in the config file or matched with a pattern
        if 'skip' in config:
            matched_pattern = input in config['skip']
            if not matched_pattern:
                for skipped_files in config['skip']:
                    # check the input script name against a pattern, e.g. in.*_imd*
                    if ('*' in skipped_files) and fnmatch.fnmatch(input, skipped_files):
                        matched_pattern = True
                        break

            if matched_pattern == True:
                msg = "   + " + input + f" ({test_id+1}/{num_tests}): skipped as specified in {configFileName}"
                print(msg)
                logger.info(msg)
                result.status = 'skipped as specified in the test configuration'
                record(result)
                test_id = test_id + 1
                continue

        # also skip if the test already completed as marked in the progress file
        if input in last_progress:
            status = str(last_progress[input].get('status', ''))
            if 'completed' in status or 'numerical checks skipped' in status:
                msg = "  + " + input + f" ({test_id+1}/{num_tests}): marked as completed or numerical checks skipped (see {progress_file})"
                logger.info(msg)
                print(msg)
                # No need to write to progress again that the run is completed
                result.status = 'skipped, already completed in a previous run'
                record(result, write_progress=False)
                test_id = test_id + 1
                continue

            if 'package not installed' in status:
                msg = "  + " + input + f" ({test_id+1}/{num_tests}): due to package not installed (see {progress_file})"
                logger.info(msg)
                print(msg)
                # No need to write to progress again that the run gets error due to missing packages
                result.status = 'skipped, package not installed in a previous run'
                record(result, write_progress=False)
                test_id = test_id + 1
                continue

        str_t = "   + " + input_test + f" ({test_id+1}/{num_tests})"
        logger.info(str_t)
        print(str_t)

        # check if a reference log file exists in the current folder: log.DDMMMYY.basename.g++.[nprocs]
        # assuming that input file names start with "in." (except in.disp, in.disp2 and in.dos in phonon/)
        basename = input_test[3:]
        ref_logfile_exist = False
        thermo_ref_file = ''

        # if there are multiple log files for different number of procs, pick the maximum number
        max_np = 1
        for file in sorted(glob.glob("log.*")):
            # looks for pattern log.{date}.{basename}.{log_compiler}.{nprocs}: log.[date].min.box.[compiler]].* vs log.[date].min.[compiler].*
            # skip over log files from previous test runs (log.{basename}.{nprocs}) and log.lammps
            parts = file.split('.')
            if len(parts) < 5:
                continue
            # get the date from the log files
            date = parts[1]
            log_compiler = file.rsplit('.',2)[1]
            pattern = f'log.{date}.{basename}.{log_compiler}.*'
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

        # if the nprocs value in the configuration file is empty then use max_np for this particular input script
        if config['nprocs'] == "":
            config['nprocs'] = str(max_np)
        else:
            # otherwise use the nprocs value in the configuration file (4 for most examples)
            logger.info(f"     Using {config['nprocs']} nprocs for {input_test} as enforced in the config file.")
            logger.info(f"     WARNING: The maximum number of procs found from the reference log files is {max_np}.")

        # store the value of nprocs to name the generated log file
        nprocs = int(config['nprocs'])

        # if valgrind is used for mem check, the run command will be
        #   mpirun -np 1 valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes /path/to/lmp_binary -in in.txt
        # so both mpiexec_numproc_flag and nprocs are empty
        if use_valgrind == True:
            config['nprocs'] = ""
            config['mpiexec_numproc_flag'] = ""
            nprocs = 1

        # walltime =   -2: skipped tests
        #              -1: failed tests
        #            >= 0: walltime in seconds (e.g. in.melt walltime = 0.2 seconds)
        # default walltime value of failed tests
        result.time = -1.0

        # run the LAMMPS binary with the input script
        status = execute(lmp_binary, config, input_test)
        output = status['stdout']
        error = status['stderr']
        returncode = int(status['returncode'])
        logfilename = status['logfilename']
        result.timeout = status['timedout']

        # restore the nprocs value in the configuration
        config['nprocs'] = saved_nprocs

        # check if the output contains ERROR
        #   there might not be a log file generated at this point, or the log file contains only the date line
        if "ERROR" in output:
            error_line = ""
            for line in output.split('\n'):
                if "ERROR" in line:
                    error_line = line
                    break
            logger.info(f"     The run terminated with {input_test} gives the following output:")
            logger.info(f"       {error_line}")
            if "Unrecognized" in output:
                result.status = f"failed, unrecognized command, package not installed, {shorten(error_line)}"
            elif "Unknown" in output:
                result.status = f"failed, unknown command, package not installed, {shorten(error_line)}"
            else:
                result.status = f"failed, {shorten(error_line)}"
            result.category = 'error'
            result.message = error_line

            logger.info(f"     Output:")
            logger.info(f"     {output}")
            logger.info(f"     Failed with {input_test}.\n")
            print(f"{result.status}")

            record(result)
            test_id = test_id + 1
            continue

        # check if a log file log.{basename}.{nprocs} exists in the current folder
        if os.path.isfile(logfilename) == False:
            msg = f"    failed, no {logfilename} generated with {input_test} with return code {returncode}.\n"
            print(msg)
            logger.info(msg)
            logger.info(f"    Output:")
            logger.info(f"    {output}")
            logger.info(f"    Error:\n{error}")

            result.category = 'error'
            if result.timeout:
                result.status = f"failed, no log file generated, {shorten(error)}"
            else:
                result.status = f"failed, no log file generated with return code {returncode}"
            result.message = shorten(error, maxlen=1000)
            record(result)
            test_id = test_id + 1
            continue
        else:
            # generate a new log file whose name has the format of log.{date}.{basename}.{compiler}.{nprocs}
            if genref == True:
                dmy = datetime.datetime.now()
                date = dmy.strftime("%d%b%y")
                shutil.copy(f"log.{basename}.{nprocs}", f"log.{date}.{basename}.{compiler}.{nprocs}")

        # parse the total wall time from the screen output (may be absent if the run did not complete)
        # NOTE: Total wall time could be 00:00:00 whereas Loop time is non-zero seconds
        walltime = None
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
                result.time = walltime
                break

        # if skip numerical checks, then skip the rest
        if skip_numerical_check == True:
            msg = "completed, skipping numerical checks"
            result.category = 'completed'
            if use_valgrind == True:
                if "All heap blocks were freed" in error:
                    msg += ", no memory leak"
                else:
                    msg += ", memory leaks detected"
                    result.category = 'failed'
                    result.memleak = True
            result.status = msg
            record(result, walltime_norm=walltime_norm)
            test_id = test_id + 1
            continue

        # if there is no ERROR in the output, but there is no Total wall time printed out
        if walltime is None:
            msg = f"     failed, no Total wall time in the output.\n"
            print(msg)
            logger.info(msg)
            logger.info(f"\n{input_test}:")
            logger.info(f"\n    Output:\n{output}")
            logger.info(f"\n    Error:\n{error}")

            result.category = 'error'
            result.status = f"failed, no Total wall time in the output, {shorten(error)}"
            result.message = shorten(error, maxlen=1000)
            record(result)
            test_id = test_id + 1
            continue

        # if there is no Step or no Loop printed out
        if "Step" not in output or "Loop" not in output:
            msg = f"    completed, but no Step nor Loop in the output.\n"
            print(msg)
            logger.info(msg)
            logger.info(f"\n{input_test}:")
            logger.info(f"\n    Output:\n{output}")
            logger.info(f"\n    Error:\n{error}")

            result.category = 'completed'
            result.status = 'completed, but no Step nor Loop in the output'
            record(result, walltime_norm=walltime_norm)
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
            result.category = 'completed'
            if use_valgrind == True:
                if "All heap blocks were freed" in error:
                    msg += ", no memory leak"
                else:
                    msg += ", memory leaks detected"
                    result.category = 'failed'
                    result.memleak = True

            result.status = msg + f", error parsing {logfilename} into YAML"
            if verbose == True:
                print(result.status)

            record(result, walltime_norm=walltime_norm)
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
                # the thermo_ref dictionary is empty
                logger.info(f"    failed, error parsing the reference log file {thermo_ref_file}.")
                result.category = 'completed'
                result.status = 'completed, numerical checks skipped, unsupported log file format'
                record(result, walltime_norm=walltime_norm)
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
                result.category = 'completed'
                result.status = 'completed, numerical checks skipped due to missing the reference log file'
                record(result, walltime_norm=walltime_norm)
                test_id = test_id + 1
                continue

        logger.info(f"     Comparing thermo output from {logfilename} against the reference log file {thermo_ref_file}")

        # check if the number of runs matches with that in the reference log file
        # maybe due to some changes to the input where the ref log file is not updated yet
        if num_runs != num_runs_ref:
            logger.info(f"     ERROR: Number of runs in {logfilename} ({num_runs}) is different from that in the reference log ({num_runs_ref})."
                        " Check README in the folder, possibly due to using mpirun with partitions or parsing the wrong reference log file.")
            result.category = 'failed'
            result.status = "failed, incomplete runs"
            record(result, walltime_norm=walltime_norm)
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
            result.category = 'failed'
            result.status = "failed, mismatched columns in the log files"
            record(result, walltime_norm=walltime_norm)
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
        mismatched_num_steps = False

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
            nthermo_steps_ref = len(thermo_ref[irun]['data'])

            if nthermo_steps_ref != nthermo_steps:
                logger.info(f"     failed: Number of thermo steps in {logfilename} ({nthermo_steps})")
                logger.info(f"     is different from that in the reference log ({nthermo_steps_ref}) in run {irun}.")
                mismatched_num_steps = True   
                continue

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
            msg = f"     mismatched columns in the log files after the first run. Check both log files for more details."
            print(msg)
            logger.info(msg)
            result.category = 'failed'
            result.status = "failed, thermo checks skipped due to mismatched log files after the first run"
            record(result, walltime_norm=walltime_norm)
            test_id = test_id + 1
            continue

        # some runs that involve the minimize command that leads to different number of steps vs the reference log file
        if mismatched_num_steps == True:
            msg = f"     mismatched num steps in the log files. Check both log files for more details."
            print(msg)
            logger.info(msg)
            result.category = 'failed'
            result.status = "failed, thermo checks skipped due to mismatched number of steps in the log files"
            record(result, walltime_norm=walltime_norm)
            test_id = test_id + 1
            continue

        if num_abs_failed > 0:
            msg = f"     {num_abs_failed} abs diff checks failed."
            print(msg)
            logger.info(msg)
            for out in failed_abs_output:
                logger.info(f"        - {out}")

            if verbose == True:
                for out in failed_abs_output:
                    print(f"        - {out}")

        if num_rel_failed > 0:
            msg = f"     {num_rel_failed} rel diff checks failed."
            print(msg)
            logger.info(msg)
            for out in failed_rel_output:
                logger.info(f"        - {out}")

            if verbose == True:
                for out in failed_rel_output:
                    print(f"        - {out}")

        result.checks = num_checks
        if num_abs_failed == 0 and num_rel_failed == 0:
            msg = f"     all {num_checks} checks passed."
            print(msg)
            logger.info(msg)
            result.category = 'passed'
            msg = "completed"
        else:
            result.category = 'failed'
            result.message = '\n'.join(failed_abs_output + failed_rel_output)
            msg = f"completed, {num_abs_failed} abs diff and {num_rel_failed} rel diff checks failed"

        # check if memleak detects from valgrind run (need to replace "mpirun" -> valgrind --leak-check=yes mpirun")
        if use_valgrind == True:
            if "All heap blocks were freed" in error:
                msg += ", no memory leak"
            else:
                msg += ", memory leaks detected"
                result.category = 'failed'
                result.memleak = True
        result.status = msg

        failed_checks = { 'num_checks': num_checks,
                          'abs_diff_failed': num_abs_failed,
                          'rel_diff_failed': num_rel_failed }
        record(result, walltime_norm=walltime_norm, failed_checks=failed_checks)
        test_id = test_id + 1

    # collect statistics over the results recorded in this call
    stat = { 'num_completed': 0,
             'num_passed': 0,
             'num_skipped': 0,
             'num_error': 0,
             'num_timeout': 0,
             'num_failed': 0,
             'num_memleak': 0,
           }
    for result in results[num_results_initial:]:
        if result.category in ('passed', 'failed', 'completed'):
            stat['num_completed'] += 1
        if result.category == 'passed':
            stat['num_passed'] += 1
        elif result.category == 'failed':
            stat['num_failed'] += 1
        elif result.category == 'error':
            stat['num_error'] += 1
        elif result.category == 'skipped':
            stat['num_skipped'] += 1
        if result.timeout:
            stat['num_timeout'] += 1
        if result.memleak:
            stat['num_memleak'] += 1
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
       - returncode:  error code returned by the process
       - timedout:    whether the run was killed after exceeding the timeout
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

    # launch the run in its own process group, so that on a timeout the whole group
    # (the shell, mpirun, and the LAMMPS processes) can be killed; with subprocess.run()
    # only the shell would be killed, leaving orphaned LAMMPS processes behind
    p = subprocess.Popen(cmd_str, shell=True, text=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, start_new_session=True)
    try:
        stdout, stderr = p.communicate(timeout=timeout)
        status = {
            'cmd_str': cmd_str,
            'stdout': stdout,
            'stderr': stderr,
            'returncode': p.returncode,
            'timedout': False,
            'logfilename': logfilename,
        }
        return status

    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(p.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            p.kill()
        p.communicate()
        msg = f"     Timeout for: {cmd_str} ({timeout}s expired)"
        logger.info(msg)
        print(msg)

    error_str = f"timeout ({timeout}s expired)"
    status = {
        'cmd_str': cmd_str,
        'stdout': "",
        'stderr': error_str,
        'returncode': -1,
        'timedout': True,
        'logfilename': logfilename,
    }
    return status

'''
   get the reference walltime by running the lmp_binary with config with an input script in the bench/ folder
      in.lj is suitable as it doesn't need any potential file, nor any extra packages
'''
def get_reference_walltime(lmp_binary, config, example_toplevel):
    cmd_str = ""
    # check if mpiexec/mpirun is used
    if config['mpiexec']:
        cmd_str += config['mpiexec'] + " " + config['mpiexec_numproc_flag'] + " " + config['nprocs'] + " "

    # guess the build folder path
    lmp_build_folder = lmp_binary.rsplit('/', 1)[0]

    # guess the bench folder
    lmp_bench_folder = example_toplevel + "/../bench/"

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

    looptime = 1.0
    for line in output.split('\n'):
        if "Total wall time" in line:
            walltime_str = line.split('time:')[1]
            hms = walltime_str.split(':')
            hours = float(hms[0])
            minutes = float(hms[1])
            seconds = float(hms[2])
            walltime = hours * 3600.0 + minutes * 60.0 + seconds
        if "Loop time" in line:
            looptime_str = line.split(' ')[3]
            seconds = float(looptime_str)
            looptime = seconds

    # there is case where total walltime with in.lj is reported as zero seconds, then use loop time
    if float(walltime) < float(config['epsilon']):
        walltime = looptime

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
    parser.add_argument("--quick-max", dest="quick_max", default=0,
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
                        num_input = len(glob.glob(f"{folder}/in.*"))
                        f.write(f"{folder} {num_input}\n")
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

            # allow to select randomly some input scripts at this point if quick_max is set
            if quick_max > 0 and len(example_inputs) > quick_max:
                example_inputs = random.sample(example_inputs, quick_max)
                msg = "\nTesting " + str(quick_max) + " randomly selected inputs"
                print(msg)
                logger.info(msg)

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
    pwd = os.getcwd()
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
            if last_progress is None:
                last_progress = {}
        except Exception:
            print(f"    Cannot open progress file {progress_file_abs} to resume, rerun all the tests")

    # get a reference walltime
    walltime_ref = get_reference_walltime(lmp_binary, config, example_toplevel)

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
    timeout_tests = 0
    failed_tests = 0
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

            if os.path.exists(directory) is False:
                continue

            # change to the directory where the input script and data files are located
            print("-"*80)
            print("Entering " + directory)
            logger.info("Entering " + directory)
            os.chdir(directory)

            all_input_list = sorted(glob.glob("in.*"))

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
            timeout_tests += stat['num_timeout']
            failed_tests += stat['num_failed']
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
        stat = iterate(lmp_binary, pwd, input_list, config,
                       results, progress_file_abs, failure_file_abs, walltime_ref, verbose, last_progress)

        completed_tests = stat['num_completed']
        skipped_tests = stat['num_skipped']
        passed_tests = stat['num_passed']
        error_tests = stat['num_error']
        timeout_tests += stat['num_timeout']
        failed_tests = stat['num_failed']
        memleak_tests = stat['num_memleak']

        all_results.extend(results)

    # print out summary:
    #  error_tests = number of runs that errored out
    #  failed_tests = number of runs that failed the numerical checks, including missing the reference log files, different num runs and num steps in a run
    #  completed_tests = number of runs that reached the end (Total wall time printed out) = failed_sests + passed_tests

    msg = "\nSummary:\n"
    msg += f"  Total number of input scripts: {total_tests}\n"
    msg += f"  - Skipped  : {skipped_tests}\n"
    msg += f"  - Error    : {error_tests}\n"
    msg += f"     - timeout  : {timeout_tests}\n"
    msg += f"  - Completed: {completed_tests}\n"
    msg += f"     - failed   : {failed_tests}\n"

    # print notice to GitHub
    if 'GITHUB_STEP_SUMMARY' in os.environ:
        with open(os.environ.get('GITHUB_STEP_SUMMARY'), 'a') as f:
            print(f"Total: {total_tests}  Skipped: {skipped_tests}  Error: {error_tests}  Timeout: {timeout_tests}"
                  f"  Failed: {failed_tests}  Passed: {passed_tests}  Completed: {completed_tests}", file=f)

    if 'valgrind' in config['mpiexec']:
        msg += f"     - memory leak detected  : {memleak_tests}\n"
    if passed_tests <= completed_tests:
        msg += f"     - numerical tests passed: {passed_tests}\n"
    msg += "\nOutput:\n"
    msg += f"  - List of failed inputs         : {failure_file}\n"
    msg += f"  - Status of the tested inputs   : {progress_file}\n"
    msg += f"  - Running log with screen output: {log_file}\n"
    msg += f"  - Testing result in JUnit XML   : {output_file}\n"

    print(msg)

    # generate a JUnit XML file with all test results for downstream reporting
    properties = {
        'operating_system': operating_system,
        'git_info': GitInfo,
        'compiler': compiler_full,
        'compile_flags': compile_flags,
        'config_file': os.path.basename(configFileName),
        'lmp_binary': lmp_binary,
    }
    write_junit_xml(output_file, all_results, suite_name=os.path.basename(configFileName),
                    properties=properties)
