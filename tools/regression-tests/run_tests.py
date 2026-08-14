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
    + distributing the input list across multiple processes by splitting the list of input scripts
      into separate runs (there are ~800 input scripts under the top-level examples) taking the
      estimated cost of the individual runs into account, so that all workers take about the same time
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
                     with error messages (for failed runs), the walltime reported by LAMMPS and the
                     measured wall-clock time of the run (both in seconds)
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
                --analyze --num-workers=8 --cost-file=costs.json

       The output of this run is 8 files folder-list-[0-7].txt that lists the subfolders
       and 8 files input-list-[0-7].txt that lists the input scripts under the top-level example folders.
       With these lists, one can launch multiple instances of run_tests.py simultaneously
       each with a list of example subfolders (Case 3), or with a list of input scripts (Case 4).
       Both lists are split such that all workers are expected to take the same time, using
       the measured run times of a previous test run from the optional --cost-file, which is
       written by merge_results.py --cost-file (see the README in this folder).
'''

from argparse import ArgumentParser
import datetime
import fnmatch
import glob
import heapq
import json
import logging
import os
import random
import re
import shutil
import signal
import subprocess
import sys
import xml.etree.ElementTree as ET
from time import monotonic
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
   time    : walltime as reported by LAMMPS in seconds (-1: failed, -2: skipped)
             note that LAMMPS prints this with a resolution of one second only
   elapsed : wall-clock time in seconds spent on launching the run and waiting
             for it to finish; this is the actual cost of the test and is
             recorded also for runs that crash or time out (0.0 if not run)
   checks  : number of performed numerical checks
   diverged: number of quantities that leave their tolerance somewhere in the run
   divat   : the timestep at which the earliest of them does so (None if nothing
             deviates or if the thermo output has no Step column)
   divrow  : how many thermo outputs into the trajectory that is
   attention: description of a problem with the input script itself that a developer
             has to fix, e.g. initial velocities that depend on the number of MPI
             processes (empty if there is none)
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
    self.elapsed = 0.0
    self.checks = checks
    self.diverged = 0
    self.divat = None
    self.divrow = None
    self.attention = ""
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
        # the measured wall-clock time is the cost of the test also when the run
        # crashed or timed out, while the walltime reported by LAMMPS is only
        # available for completed runs and has a resolution of one second
        case.set('time', f"{max(float(result.elapsed), 0.0):.3f}")
        total_time += max(float(result.elapsed), 0.0)

        # where in the run the output starts to deviate from the reference log file;
        # this tells apart a run that is wrong from the start from one that slowly
        # drifts apart, which need very different responses
        if result.diverged:
            case.set('diverged', str(result.diverged))
            case.set('diverged-row', str(result.divrow))
            if result.divat is not None:
                case.set('diverged-at', str(result.divat))
        if result.attention:
            case.set('attention', shorten(result.attention))

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
        elif result.message:
            # a test that passed, but has something worth reporting, e.g. quantities
            # that deviate from the reference before the last thermo output
            ET.SubElement(case, 'system-out').text = result.message

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

    # reference log files in this folder that belong to none of its input scripts, so
    # that an input without a reference log file can point at the likely reason
    orphan_logs = orphaned_reference_logs('.')
    if orphan_logs:
        logger.info(f"     {len(orphan_logs)} reference log file(s) in {input_folder} belong to"
                    f" no input script: {' '.join(orphan_logs)}")

    EPSILON = np.float64(config['epsilon'])
    nugget = float(config['nugget'])
    genref = config['genref']
    compiler = config['compiler']
    use_valgrind = 'valgrind' in config['mpiexec']

    # inputs without a reference log file are run with their run commands shortened to
    # this many steps, since all that can be checked for them is that they do not crash
    smoke_steps = int(config['smoke_steps']) if str(config.get('smoke_steps', "")) != "" else 0

    # runs are killed after this many seconds (must match the default in execute());
    # an input that needs a good fraction of it is a production run rather than an
    # example and should be shortened in the repository
    timeout_seconds = int(config['timeout']) if str(config.get('timeout', "")) != "" else 60
    long_run_seconds = 0.25 * timeout_seconds

    # record the outcome of a test: append it to the results list, write one line
    # to the progress file and, for failed or errored tests, to the failure file;
    # each line is a flow-style YAML mapping so that the whole file parses as a
    # single YAML mapping keyed by the input script names
    def record(result, walltime_norm=None, failed_checks=None, write_progress=True):
        results.append(result)
        if not write_progress:
            return
        entry = { 'folder': result.folder, 'status': result.status }
        if result.attention:
            entry['attention'] = result.attention
        if failed_checks is not None:
            entry['failed_checks'] = failed_checks
        entry['walltime'] = float(result.time)
        if walltime_norm is not None:
            entry['walltime_norm'] = float(walltime_norm)
        entry['elapsed'] = round(float(result.elapsed), 3)
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

        # input scripts that need a multi-partition run write one log file per partition
        # and cannot be tested with a configuration that runs them in a single partition
        if not ({'-p', '-partition'} & set(config['args'].split())):
            command = needs_partitions(input)
            if command:
                msg = ("   + " + input + f" ({test_id+1}/{num_tests}): skipped, needs a"
                       f" multi-partition run because of \"{command}\"")
                print(msg)
                logger.info(msg)
                result.status = f'skipped, needs a multi-partition run ("{command}")'
                record(result)
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
        for nprocs_ref, file in find_reference_logs('.', basename):
            # if using valgrind or running in serial, then use the log file with 1 proc
            if use_valgrind == True or config['mpiexec'] == "":
                if nprocs_ref == 1:
                    max_np = nprocs_ref
                    ref_logfile_exist = True
                    thermo_ref_file = file
                    break
            else:
                if max_np <= nprocs_ref:
                    max_np = nprocs_ref
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

        # The initial velocities of some inputs depend on the domain decomposition, so
        # their output can never match a reference log file that was written with a
        # different number of MPI processes.  This has to be corrected in the input
        # script, after which its reference log files have to be recreated.  The number
        # of MPI processes of the reference log file is read from the log file itself:
        # the number in its name counts the OpenMP threads as well.
        if ref_logfile_exist:
            ref_nprocs, ref_nthreads = reference_log_decomposition(thermo_ref_file)
            if ref_nprocs is not None and ref_nprocs != nprocs:
                reason = velocity_decomposition_dependent(input_test)
                if reason:
                    result.attention = (f"{reason}: cannot match {thermo_ref_file}, which was"
                                        f" written with {ref_nprocs} MPI processes instead of"
                                        f" {nprocs}")
                    logger.info(f"     {result.attention}")

        # walltime =   -2: skipped tests
        #              -1: failed tests
        #            >= 0: walltime in seconds (e.g. in.melt walltime = 0.2 seconds)
        # default walltime value of failed tests
        result.time = -1.0

        # Without a reference log file nothing can be checked, so the only thing left to
        # test is that the input starts and runs without crashing.  There is no reason to
        # run it to the end for that, so a copy with shortened run commands is used.  This
        # is where most of the time of a full test run used to go: those are the inputs
        # that nobody trimmed for testing and that hit the timeout.
        smoke_input = ""
        shortened = 0
        smoke_failed = False
        if (smoke_steps > 0) and (ref_logfile_exist == False) and (use_valgrind == False) \
                and (genref == False) and not os.path.isfile('thermo.' + input + '.yaml'):
            smoke_input = input_test + '.smoke'
            try:
                shortened = write_shortened_input(input_test, smoke_input, smoke_steps)
            except OSError as err:
                logger.info(f"     Cannot write {smoke_input}: {err}")
                smoke_input = ""
            if smoke_input and shortened == 0:
                # nothing to shorten, so run the input script itself
                try:
                    os.remove(smoke_input)
                except OSError:
                    pass
                smoke_input = ""
            elif smoke_input:
                logger.info(f"     No reference log file: running {input_test} with"
                            f" {shortened} run command(s) shortened to {smoke_steps} steps")

        # run the LAMMPS binary with the input script
        status = execute(lmp_binary, config, smoke_input if smoke_input else input_test)

        # Shortening the runs can break an input script, for example when a variable
        # refers to a fix that only produces output after more steps than the shortened
        # run has.  There is no way to tell such an artifact from a real problem, so the
        # input script is run unchanged in that case, as it was before.  A shortened run
        # that hits the timeout is not such a case: its cost is in the setup rather than
        # in the run commands, and running it unchanged would only waste the timeout twice.
        if smoke_input and not status['timedout'] \
                and (("ERROR" in status['stdout']) or (status['returncode'] != 0)):
            logger.info(f"     {input_test} does not work with shortened runs, running it unchanged")
            elapsed = status['elapsed']
            try:
                os.remove(smoke_input)
            except OSError:
                pass
            smoke_input = ""
            smoke_failed = True
            status = execute(lmp_binary, config, input_test)
            status['elapsed'] += elapsed

        output = status['stdout']
        error = status['stderr']
        returncode = int(status['returncode'])
        logfilename = status['logfilename']
        result.timeout = status['timedout']
        result.elapsed = status['elapsed']

        # Many example inputs are older contributions that were never trimmed for
        # testing: they run a production number of steps, which costs most of the time
        # of a full test run and, for the ones that hit the timeout, produces no result
        # at all.  Those have to be shortened in the repository, so flag them.
        if result.timeout or result.elapsed >= long_run_seconds:
            note = (f"the run needs {result.elapsed:.0f} s"
                    + (f" and hits the timeout of {timeout_seconds} s" if result.timeout else "")
                    + ": a production sized run that should be shortened in the repository")
            if smoke_failed:
                note += " (shortening its run commands makes it fail, so it has to be done by hand)"
            elif smoke_input:
                note += (f" (it does not even finish with its run commands shortened to"
                         f" {smoke_steps} steps, so the cost is in the setup)")
            result.attention = f"{result.attention}; {note}" if result.attention else note
            logger.info(f"     {note}")

        # restore the nprocs value in the configuration
        config['nprocs'] = saved_nprocs

        # the shortened copy of the input script is not needed any longer
        if smoke_input:
            try:
                os.remove(smoke_input)
            except OSError:
                pass

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
            # utils::check_packages_for_style() names the package a style belongs to, so
            # "Unrecognized ... style" without that part means that the style does not
            # exist in LAMMPS at all and the input script uses an outdated name
            if "package which is not enabled" in error_line:
                result.status = f"failed, package not installed, {shorten(error_line)}"
            elif "missing because of a dependency" in error_line:
                result.status = f"failed, style missing from this build, {shorten(error_line)}"
            elif "Unrecognized" in error_line and " style " in error_line:
                result.status = f"failed, no such style in this LAMMPS version, {shorten(error_line)}"
                result.attention = ("names a style that does not exist in any package, so this is"
                                    " a typo or an outdated name that has to be fixed in the input")
            elif "Unknown command" in error_line:
                # a command of a package that is not installed, or not a LAMMPS input at all
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
                if smoke_input:
                    result.status = (f"completed, no reference log file, only checked that a run"
                                     f" shortened to {smoke_steps} steps does not crash")
                else:
                    result.status = 'completed, numerical checks skipped due to missing the reference log file'
                # a reference log file in this folder that belongs to no input script is
                # the likely reason: one of the two names has a typo or was not renamed
                if orphan_logs:
                    note = ("no reference log file matches this input, while "
                            + ", ".join(orphan_logs[:3])
                            + (", ..." if len(orphan_logs) > 3 else "")
                            + " in this folder match no input script at all: check the names")
                    result.attention = f"{result.attention}; {note}" if result.attention else note
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

        # check if overrides for this input scipt is specified
        overrides = {}
        if 'overrides' in config:
            if input_test in config['overrides']:
                overrides = config['overrides'][input_test]

        # compare the output against the reference values, one row of thermo output
        # at a time, and decide the verdict from the last row of each run
        comparison = compare_thermo(thermo, thermo_ref, config['tolerance'], overrides,
                                    EPSILON, nugget)

        # one of the log files has no usable thermo output for this run
        if comparison['status'] == 'no thermo data':
            msg = (f"     no thermo data to compare in run {comparison['run']} of"
                   f" {logfilename} or {thermo_ref_file}.")
            print(msg)
            logger.info(msg)
            result.category = 'completed'
            result.status = 'completed, numerical checks skipped, unsupported log file format'
            record(result, walltime_norm=walltime_norm)
            test_id = test_id + 1
            continue

        # the log files do not have the same number of thermo columns
        if comparison['status'] == 'mismatched columns':
            msg = (f"     mismatched columns in the log files in run {comparison['run']}."
                   " Check both log files for more details.")
            print(msg)
            logger.info(msg)
            result.category = 'failed'
            result.status = "failed, thermo checks skipped due to mismatched log files after the first run"
            record(result, walltime_norm=walltime_norm)
            test_id = test_id + 1
            continue

        # some runs that involve the minimize command that leads to different number of steps vs the reference log file
        if comparison['status'] == 'mismatched steps':
            msg = (f"     mismatched num steps in the log files in run {comparison['run']}."
                   " Check both log files for more details.")
            print(msg)
            logger.info(msg)
            result.category = 'failed'
            result.status = "failed, thermo checks skipped due to mismatched number of steps in the log files"
            record(result, walltime_norm=walltime_norm)
            test_id = test_id + 1
            continue

        num_checks = comparison['num_checks']
        num_abs_failed = comparison['num_abs_failed']
        num_rel_failed = comparison['num_rel_failed']
        failed_abs_output = comparison['failed_abs']
        failed_rel_output = comparison['failed_rel']

        # judge the earliest deviation: too early for the chaotic nature of a classical
        # MD trajectory to explain it, or late enough that it is expected
        assessment = divergence_assessment(comparison['diverged_at'], comparison['diverged_row'],
                                           comparison['thermo_every'])
        if assessment and result.attention:
            assessment += f"; NOTE: {result.attention}"
        divergence = ([assessment] if assessment else []) + format_divergence(comparison['diverged'])

        if verbose == True:
            width = 20
            print("        Quantities".ljust(width) + "Output".center(width) + "Reference".center(width) +
                  "Abs Diff Check".center(width) + "Rel Diff Check".center(width))
            for entry in comparison['last_row']:
                abs_check = "PASSED" if entry['abs_ok'] else "FAILED"
                rel_check = "PASSED" if entry['rel_ok'] else "FAILED"
                print(f"        {entry['quantity'].ljust(width)} {str(entry['value']).rjust(20)}"
                      f" {str(entry['reference']).rjust(20)} {abs_check.rjust(20)} {rel_check.rjust(20)}")

        if comparison['skipped_rows']:
            logger.info(f"     skipped {comparison['skipped_rows']} lines that look like thermo"
                        f" output but do not consist of numbers in {logfilename} or {thermo_ref_file}")

        if divergence:
            logger.info(f"     {len(comparison['diverged'])} quantities deviate from the reference log file"
                        f" in {comparison['num_rows']} compared thermo outputs:")
            for line in divergence:
                logger.info(f"        - {line}")
            if verbose == True:
                for line in divergence:
                    print(f"        - {line}")

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

        # quantities that deviate from the reference somewhere in the run, but are back
        # within their tolerance at the last thermo output where the verdict is decided
        transient = [entry for entry in comparison['diverged']
                     if entry['last_abs'] <= entry['abs_tol'] and entry['last_rel'] <= entry['rel_tol']]

        result.checks = num_checks
        if num_abs_failed == 0 and num_rel_failed == 0:
            msg = f"     all {num_checks} checks passed."
            print(msg)
            logger.info(msg)
            result.category = 'passed'
            msg = "completed"
            if transient:
                msg += (f", but {len(transient)} quantities deviate from the reference"
                        " before the last thermo output")
                result.message = '\n'.join(divergence)
        else:
            result.category = 'failed'
            # the divergence table comes first: it is what tells apart a real
            # difference from a run that only drifts apart towards its end, and the
            # report may be truncated by the downstream tools
            result.message = '\n'.join(divergence + [''] + failed_abs_output + failed_rel_output)
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
        if comparison['diverged']:
            # the earliest deviation from the reference log file: a deviation in the
            # first thermo outputs is a difference in what is computed, while one after
            # a few thousand steps is the chaotic nature of a classical MD trajectory
            result.diverged = len(comparison['diverged'])
            result.divat = comparison['diverged_at']
            result.divrow = comparison['diverged_row']
            failed_checks['diverged'] = result.diverged
            failed_checks['diverged_row'] = result.divrow
            if result.divat is not None:
                failed_checks['diverged_at'] = result.divat
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

# commands that only work when LAMMPS is started with more than one partition, i.e. with
# "-partition NxM" and a matching number of MPI processes.  Without that they either stop
# with an error or compute something meaningless, and their thermo output goes to one log
# file per partition instead of the single log file that is compared here.  The list
# follows the universe->nworlds checks in the sources.
PARTITION_COMMANDS = [
    re.compile(r'^\s*(neb|neb/spin|prd|tad|temper|temper/npt|temper/grem)\s', re.IGNORECASE),
    re.compile(r'^\s*fix\s+\S+\s+\S+\s+(pimd/langevin|pimd/nvt|alchemy|gemc|neb|neb/spin)\s',
               re.IGNORECASE),
    re.compile(r'^\s*run_style\s+verlet/split\s*$', re.IGNORECASE),
    re.compile(r'^\s*kspace_style\s+pppm/rk\s', re.IGNORECASE),
    re.compile(r'^\s*partition\s', re.IGNORECASE),
    re.compile(r'^\s*variable\s+\S+\s+(universe|uloop)\s', re.IGNORECASE),
]

'''
    check whether an input script needs a multi-partition run

    return the offending command, or an empty string if the input script can be run
    with a single partition
'''
def needs_partitions(input_file):
    try:
        with open(input_file, 'r', errors='ignore') as f:
            for line in f:
                if line.lstrip().startswith('#'):
                    continue
                for pattern in PARTITION_COMMANDS:
                    match = pattern.match(line)
                    if match:
                        return ' '.join(match.group(0).split())
    except OSError:
        pass
    return ""

'''
    find the reference log files in a folder that belong to none of its input scripts

    A reference log file is only found when its name follows log.{date}.{basename}.
    {compiler}.{nprocs} with the basename of an input script in the same folder.  A
    typo in either name, or renaming one without the other, silently turns the input
    into one that cannot be checked against anything.

    return the list of those file names, without the log files of multi-partition runs
'''
def orphaned_reference_logs(folder):
    claimed = set()
    for path in sorted(glob.glob(os.path.join(folder, 'in.*'))):
        for nprocs, log in find_reference_logs(folder, os.path.basename(path)[3:]):
            claimed.add(log)

    orphans = []
    for path in sorted(glob.glob(os.path.join(folder, 'log.*'))):
        log = os.path.basename(path)
        parts = log.split('.')
        if len(parts) < 5 or not parts[-1].isdigit() or log in claimed:
            continue
        # log.{date}.{basename}.{compiler}.{nprocs}.{partition} of a multi-partition run
        if parts[-2].isdigit():
            continue
        orphans.append(log)
    return orphans

'''
    read the domain decomposition from the header of a log file

    The number at the end of the name of a reference log file is the number of MPI
    processes times the number of OpenMP threads per process, which is what LAMMPS
    reports in its "Loop time of ... on N procs" line.  A run with 2 MPI processes and
    2 threads each therefore has the same reference log file as one with 4 MPI
    processes, although the two decompose the simulation box differently.  What matters
    for the reproducibility of atom IDs and initial velocities is the number of MPI
    processes alone, so that is read from the log file itself.

    return a tuple of (number of MPI processes, number of threads per process), with
    either of them None if the log file does not report it
'''
def reference_log_decomposition(filename):
    grid = re.compile(r'^\s*(\d+) by (\d+) by (\d+) MPI processor grid')
    omp = re.compile(r'^\s*using (\d+) OpenMP thread')
    nprocs = None
    nthreads = None
    try:
        with open(filename, 'r', errors='ignore') as f:
            for line in f:
                match = grid.match(line)
                if match and nprocs is None:
                    nprocs = int(match.group(1)) * int(match.group(2)) * int(match.group(3))
                match = omp.match(line)
                if match and nthreads is None:
                    nthreads = int(match.group(1))
                # both are printed in the header, no need to read the whole log file
                if nprocs is not None and nthreads is not None:
                    break
    except OSError:
        pass
    return nprocs, nthreads

'''
    check whether the initial velocities of an input script depend on the number of
    MPI processes, which makes its output impossible to compare against a reference log
    file that was created with a different number of processes

    Following the documentation of the velocity command:
      + "loop geom" seeds the random numbers with the coordinates of each atom and always
        assigns the same velocity to an atom
      + "loop local" adjusts the seed per process and always gives different velocities
      + "loop all" (the default) assigns velocities by atom ID, which only gives the same
        result when the atoms were read from a data file: atoms created with create_atoms
        get their IDs assigned depending on the domain decomposition

    Such an input needs to be corrected at the source (add "loop geom" to the velocity
    command) and its reference log files have to be recreated.

    return a description of the problem, or an empty string if the velocities are
    reproducible or if the input script does not create any
'''
def velocity_decomposition_dependent(input_file):
    create = re.compile(r'^\s*velocity\s+\S+\s+create\s', re.IGNORECASE)
    lines = []
    try:
        with open(input_file, 'r', errors='ignore') as f:
            # join the lines that are continued with a trailing "&"
            text = f.read().replace('&\n', ' ')
    except OSError:
        return ""

    lines = [line for line in text.splitlines() if not line.lstrip().startswith('#')]
    velocity = [line for line in lines if create.match(line)]
    if not velocity:
        return ""
    if all('loop geom' in ' '.join(line.split()) for line in velocity):
        return ""
    if any(re.search(r'\bloop\s+local\b', line, re.IGNORECASE) for line in velocity):
        return 'velocity create with "loop local"'
    if any(re.match(r'^\s*create_atoms\s', line, re.IGNORECASE) for line in lines):
        return 'velocity create with the default "loop all" and atoms from create_atoms'
    return ""

'''
    copy an input script and reduce the number of steps of its run commands

    This is for input scripts without a reference log file: their output cannot be
    checked against anything, so the only thing left to test is that they start and run
    without crashing, and there is no reason to run them to the end.

    input_file : the input script to read
    output_file: the shortened input script to write
    maxsteps   : upper limit for the number of steps of a run command
    return the number of run commands that were shortened
'''
def write_shortened_input(input_file, output_file, maxsteps):
    # "run 10000 upto", "run ${nsteps} post no", ...: only the number of steps is replaced
    pattern = re.compile(r'^(\s*run\s+)(\S+)(.*)$', re.IGNORECASE)
    shortened = 0
    with open(input_file, 'r', errors='ignore') as src, open(output_file, 'w') as dst:
        for line in src:
            match = pattern.match(line)
            if match:
                steps = match.group(2)
                # never make a run longer than it was, e.g. the "run 0" of a single point
                nsteps = min(int(steps), maxsteps) if steps.isdigit() else maxsteps
                if str(nsteps) != steps:
                    line = f"{match.group(1)}{nsteps}{match.group(3)}\n"
                    shortened += 1
            dst.write(line)
    return shortened

'''
    compare the thermo output of a run against the reference log file

    All thermo rows are compared, not only the last one.  Comparing only the last row
    hides where a deviation starts, which is what tells apart a run that slowly drifts
    apart from the reference from one that is wrong from the very first output.

    thermo, thermo_ref: thermo data of the run and of the reference log file
    tolerance : the "tolerance" section of the test configuration
    overrides : per-input tolerances that replace the ones in "tolerance"
    epsilon   : reference values below this are considered to be zero
    nugget    : added to the reference value when it is considered to be zero

    return a dictionary with
      status        : 'ok', 'mismatched columns', 'mismatched steps' or 'no thermo data'
      run           : the run the problem was found in (for the latter three)
      skipped_rows  : number of rows that do not consist of numbers and were skipped;
                      some log files contain lines that look like thermo output but are
                      not, e.g. the SHAKE statistics printed by fix shake
      num_checks    : number of checks performed at the last row of each run
      num_abs_failed: number of failed absolute difference checks at those rows
      num_rel_failed: number of failed relative difference checks at those rows
      failed_abs    : messages for the failed absolute difference checks
      failed_rel    : messages for the failed relative difference checks
      num_rows      : total number of compared thermo rows
      diverged      : one entry per run and quantity that exceeds its tolerance in at
                      least one row, with the first diverging row, the largest deviation
                      and the deviation in the last row
      diverged_at   : the earliest timestep at which any quantity leaves its tolerance
                      (None if the thermo output has no Step column)
      diverged_row  : how many thermo outputs into the whole trajectory that is
      thermo_every  : the number of steps between two thermo outputs of the first run,
                      which is the resolution at which a deviation can be detected at all
      last_row      : the compared values in the last row of each run, for verbose output
'''
def compare_thermo(thermo, thermo_ref, tolerance, overrides, epsilon, nugget):
    result = { 'status': 'ok', 'run': 0, 'num_checks': 0,
               'num_abs_failed': 0, 'num_rel_failed': 0,
               'failed_abs': [], 'failed_rel': [], 'num_rows': 0, 'skipped_rows': 0,
               'diverged': [], 'diverged_at': None, 'diverged_row': None,
               'thermo_every': None, 'last_row': [] }

    for irun in range(len(thermo)):
        keywords = thermo[irun]['keywords']
        data = thermo[irun]['data']
        data_ref = thermo_ref[irun]['data']
        if not keywords or not data or not data_ref:
            result['status'] = 'no thermo data'
            result['run'] = irun
            return result

        if len(keywords) != len(thermo_ref[irun]['keywords']):
            result['status'] = 'mismatched columns'
            result['run'] = irun
            return result

        if len(data) != len(data_ref):
            result['status'] = 'mismatched steps'
            result['run'] = irun
            return result

        # skip the rows that do not consist of numbers: converting a log file into
        # thermo data also picks up lines that only look like thermo output, e.g. the
        # statistics that fix shake prints
        rows = []
        for row in range(len(data)):
            try:
                [float(value) for value in data[row] + data_ref[row]]
            except (TypeError, ValueError):
                continue
            rows.append(row)
        result['skipped_rows'] += len(data) - len(rows)
        if not rows:
            result['status'] = 'no thermo data'
            result['run'] = irun
            return result

        rows_before = result['num_rows']
        result['num_rows'] += len(rows)
        last_row = rows[-1]

        # the timesteps are more useful than the row numbers when looking at the log files
        steps = None
        if 'Step' in keywords:
            steps = [row[keywords.index('Step')] for row in data]
            # how far apart two thermo outputs are: a deviation cannot be seen any
            # earlier than the first thermo output after the one where it starts
            if result['thermo_every'] is None and len(rows) > 1:
                result['thermo_every'] = steps[rows[1]] - steps[rows[0]]

        for i in range(len(keywords)):
            quantity = keywords[i]
            if quantity not in tolerance and quantity not in overrides:
                # no tolerances are defined for this quantity in the config file
                continue
            if quantity in tolerance:
                abs_tol = float(tolerance[quantity]['abs'])
                rel_tol = float(tolerance[quantity]['rel'])
            # the per-input overrides replace the global tolerance values
            if quantity in overrides:
                abs_tol = float(overrides[quantity]['abs'])
                rel_tol = float(overrides[quantity]['rel'])

            first = None
            max_rel = -1.0
            max_row = rows[0]
            for row in rows:
                val = float(data[row][i])
                ref = float(data_ref[row][i])
                abs_diff = abs(val - ref)
                if abs(ref) > epsilon:
                    rel_diff = abs_diff / abs(ref)
                else:
                    rel_diff = abs_diff / abs(ref + nugget)

                if rel_diff > max_rel:
                    max_rel = rel_diff
                    max_row = row
                if first is None and (abs_diff > abs_tol or rel_diff > rel_tol):
                    first = row

                if row == last_row:
                    # the verdict of the test is decided by the last row of the run
                    result['num_checks'] += 2
                    result['last_row'].append({
                        'quantity': quantity, 'value': val, 'reference': ref,
                        'abs_ok': abs_diff <= abs_tol, 'rel_ok': rel_diff <= rel_tol })
                    if abs_diff > abs_tol:
                        result['failed_abs'].append(
                            f"Run {irun}: {quantity}: actual ({abs_diff:0.2e}) > expected ({abs_tol:0.2e})")
                        result['num_abs_failed'] += 1
                    if rel_diff > rel_tol:
                        result['failed_rel'].append(
                            f"Run {irun}: {quantity}: actual ({rel_diff:0.2e}) > expected ({rel_tol:0.2e})")
                        result['num_rel_failed'] += 1
                    if first is not None:
                        result['diverged'].append({
                            'run': irun, 'quantity': quantity, 'rows': len(data),
                            'first_row': first, 'first_step': steps[first] if steps else None,
                            'max_rel': max_rel, 'max_row': max_row,
                            'max_step': steps[max_row] if steps else None,
                            'last_abs': abs_diff, 'last_rel': rel_diff,
                            'abs_tol': abs_tol, 'rel_tol': rel_tol,
                        })
                        # the earliest deviation over the whole trajectory: a run that
                        # continues a previous one starts out deviating if that one did,
                        # so this has to be tracked across the runs and not per run
                        row_overall = rows_before + rows.index(first)
                        if result['diverged_row'] is None or row_overall < result['diverged_row']:
                            result['diverged_row'] = row_overall
                            result['diverged_at'] = steps[first] if steps else None
    return result

# A classical MD trajectory is chaotic: the tiniest difference in the computed forces,
# as it comes from a different machine, compiler, level of optimization or number of MPI
# processes, grows exponentially.  It takes about this many steps before such a rounding
# difference has grown enough to be visible in the printed thermo output, and a few
# thousand more before it reaches the printed precision of all quantities.  Anything that
# deviates later than this is expected and says nothing about the correctness of the code.
CHAOS_STEPS = 1000

# A real difference in what is computed, on the other hand, is normally visible within
# the first steps, and in subtle cases after one or two hundred.  A deviation that shows
# up this early is worth looking at.
SUSPICIOUS_STEPS = 200

'''
    judge the earliest deviation of a run from its reference log file

    diverged_at : timestep of the earliest deviation (None if there is no Step column)
    diverged_row: how many thermo outputs into the trajectory that is
    thermo_every: number of steps between two thermo outputs

    return a short assessment for the test report, or an empty string if nothing deviates
'''
def divergence_assessment(diverged_at, diverged_row, thermo_every):
    if diverged_row is None:
        return ""
    if diverged_row == 0:
        return ("deviates in the very first thermo output, before the trajectory can"
                " diverge: a difference in the setup or in the computation")
    if diverged_at is None:
        return f"deviates after {diverged_row} thermo outputs (no Step column in the output)"

    # a deviation cannot be seen before the first thermo output after it starts, so with
    # thermo output every few thousand steps the two cases cannot be told apart
    if thermo_every and thermo_every > SUSPICIOUS_STEPS and diverged_row <= 1:
        return (f"deviates in the second thermo output (step {diverged_at}); with thermo"
                f" output only every {thermo_every} steps this cannot be told apart from"
                " chaotic divergence")
    if diverged_at <= SUSPICIOUS_STEPS:
        return (f"deviates after {diverged_at} steps, too early for chaotic divergence:"
                " a difference in the computation until proven otherwise")
    if diverged_at <= CHAOS_STEPS:
        return f"deviates after {diverged_at} steps, early for chaotic divergence"
    return f"deviates after {diverged_at} steps, consistent with chaotic divergence"

'''
    describe how a quantity deviates from the reference, to tell apart the cases that
    need a shorter run or a looser tolerance from the ones that need a closer look

    entry: one element of the "diverged" list returned by compare_thermo()
'''
def divergence_pattern(entry):
    if entry['rows'] < 2:
        return "single thermo output"
    if entry['first_row'] == 0:
        # the first thermo output of a continuation run repeats the last one of the
        # previous run, so only the very first run says something about the setup
        if entry['run'] == 0:
            return "differs already in the first thermo output, check the setup"
        return "differs already at the start of this run"
    if entry['last_rel'] >= 0.5 * entry['max_rel']:
        if entry['first_row'] > 0.5 * (entry['rows'] - 1):
            return "grows only towards the end, a shorter run may be enough"
        return "grows over the whole run"
    return "largest deviation before the end of the run"

'''
    format the divergence data from compare_thermo() as a table for the test report
'''
def format_divergence(diverged, maxlen=25):
    lines = []
    for entry in diverged[:maxlen]:
        where = (f"step {entry['first_step']}" if entry['first_step'] is not None
                 else f"row {entry['first_row']}")
        peak = (f"step {entry['max_step']}" if entry['max_step'] is not None
                else f"row {entry['max_row']}")
        lines.append(f"Run {entry['run']}: {entry['quantity']}: first deviates at {where} "
                     f"({entry['first_row'] + 1} of {entry['rows']} thermo outputs), "
                     f"largest relative deviation {entry['max_rel']:0.2e} at {peak}, "
                     f"{entry['last_rel']:0.2e} at the end "
                     f"(rel. tolerance {entry['rel_tol']:0.2e}); {divergence_pattern(entry)}")
    if len(diverged) > maxlen:
        lines.append(f"... and {len(diverged) - maxlen} more")
    return lines

'''
  find the reference log files for an input script in a folder

  The expected name of a reference log file is log.{date}.{basename}.{compiler}.{nprocs}
  where basename is the input script name without the leading "in.".  Log files written
  by previous test runs (log.{basename}.{nprocs}) and log.lammps are skipped.

  folder  : folder that contains the input script and its reference log files
  basename: input script name with the leading "in." removed
  return  : list of (nprocs, filename) tuples with the bare file names, sorted by
            file name so that callers can reproduce a deterministic selection
'''
def find_reference_logs(folder, basename):
    logs = []
    for path in sorted(glob.glob(os.path.join(folder, 'log.*'))):
        file = os.path.basename(path)
        # log.[date].min.box.[compiler].* vs log.[date].min.[compiler].*
        parts = file.split('.')
        if len(parts) < 5:
            continue
        # get the date from the log files
        date = parts[1]
        log_compiler = file.rsplit('.', 2)[1]
        pattern = f'log.{date}.{basename}.{log_compiler}.*'
        if fnmatch.fnmatch(file, pattern):
            nprocs = file.rsplit('.', 1)[1]
            if nprocs.isnumeric():
                logs.append((int(nprocs), file))
    return logs

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
       - elapsed:     wall-clock time in seconds spent on the run, including the
                      process startup and (for MPI runs) the mpirun overhead
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
    start = monotonic()
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
            'elapsed': monotonic() - start,
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
        'elapsed': monotonic() - start,
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

# COST ESTIMATES AND WORK DISTRIBUTION
#
# The cost of the individual runs spans several orders of magnitude, so splitting
# the list of input scripts into equally long chunks leaves most workers idle while
# a few of them are still running the expensive inputs.  The functions below assign
# a cost estimate to every input script and distribute the input scripts over the
# workers so that all of them are expected to finish at about the same time.

# Input scripts without any timing information get the cost of this quantile of the
# input scripts that do have one.  A pessimistic guess is better than an optimistic
# one: scheduling a short run as if it were long only shifts it to the front, while
# a long run that is scheduled last delays the whole test run.
UNKNOWN_COST_QUANTILE = 0.9

# time in seconds it takes to launch a run and to process its output; added to every
# cost estimate so that the many very short runs are not considered to be free
STARTUP_COST = 0.4

'''
    the key under which an input script is stored in the cost file: the path of the
    input script relative to the top-level examples folder, e.g. PACKAGES/tally/in.pe
    (using forward slashes, so that cost files can be shared between machines)

    the same key is used as classname/name in the JUnit XML output, which is what
    merge_results.py writes the cost file from
'''
def example_key(input_path, example_toplevel):
    path = os.path.abspath(input_path)
    if example_toplevel:
        try:
            path = os.path.relpath(path, os.path.abspath(example_toplevel))
        except ValueError:
            pass
    marker = os.sep + 'examples' + os.sep
    if marker in path:
        path = path.split(marker, 1)[1]
    return path.replace(os.sep, '/')

'''
    read a cost file as written by merge_results.py --cost-file

    returns a dictionary that maps the input script keys to the measured wall-clock
    time in seconds; a missing or unreadable file results in an empty dictionary
'''
def load_cost_file(filename):
    if not filename:
        return {}
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
    except (OSError, ValueError) as err:
        print(f"WARNING: cannot read cost file {filename}: {err}")
        return {}
    costs = data.get('costs', data) if isinstance(data, dict) else {}
    return {key: float(value) for key, value in costs.items()}

'''
    the number of MPI processes times the number of OpenMP threads per process that a
    test configuration uses, which is what the reference log files are named after and
    what the "Loop time of ... on N procs" line of a log file reports

    nprocs: number of MPI processes from the test configuration
    args  : the command line arguments from the test configuration
'''
def config_ntasks(nprocs, args):
    # "-pk omp 2" for the OPENMP package, "-k on t 2" for the KOKKOS package
    match = re.search(r'-pk\s+omp\s+(\d+)|-k\s+on\s+t\s+(\d+)', args)
    if match:
        return nprocs * int(match.group(1) or match.group(2))
    return nprocs

'''
    estimate the cost of an input script from the loop times in its reference log file

    This is only a rough hint: the reference log files were generated on different
    machines over many years.  The loop times are summed over all runs in the log file
    and scaled from the number of procs of the log file to the number of procs used for
    the test, assuming perfect scaling.

    return the estimated time in seconds, or None if there is no usable reference log
'''
def reference_log_cost(folder, basename, nprocs):
    logs = find_reference_logs(folder, basename)
    if not logs:
        return None

    # prefer the reference log file for the number of procs used for the test
    candidates = [entry for entry in logs if entry[0] == nprocs]
    if not candidates:
        candidates = [max(logs)]
    nprocs_ref, filename = candidates[-1]

    looptime = 0.0
    found = False
    try:
        with open(os.path.join(folder, filename), 'r') as f:
            for line in f:
                # "Loop time of 1.61045 on 1 procs for 250 steps with 4000 atoms"
                if line.startswith('Loop time of '):
                    try:
                        looptime += float(line.split()[3])
                        found = True
                    except (IndexError, ValueError):
                        pass
    except OSError:
        return None

    if not found:
        return None
    return looptime * float(nprocs_ref) / float(max(nprocs, 1))

'''
    assign a cost estimate to every input script in input_list

    costs           : measured wall-clock times from a previous run (see load_cost_file)
    example_toplevel: top-level examples folder the cost file keys are relative to
    nprocs          : number of procs the tests are run with
    timeout         : the runs are killed after this many seconds, so no estimate
                      taken from a reference log file may be larger than this
    return a tuple of (dictionary mapping the input script paths to their estimated
    cost in seconds, dictionary with the number of inputs per source of the estimate)
'''
def estimate_costs(input_list, example_toplevel, costs, nprocs, timeout):
    estimate = {}
    sources = { 'measured': 0, 'reference log': 0, 'unknown': 0 }
    for input_path in input_list:
        key = example_key(input_path, example_toplevel)
        if key in costs:
            estimate[input_path] = max(costs[key], 0.0)
            sources['measured'] += 1
            continue
        folder, name = os.path.split(input_path)
        cost = reference_log_cost(folder, name[3:], nprocs) if name.startswith('in.') else None
        if cost is None:
            estimate[input_path] = None
            sources['unknown'] += 1
        else:
            # some reference log files are from much longer production runs than the
            # ones in the example inputs, so the estimate is capped at the timeout
            estimate[input_path] = min(cost, float(timeout))
            sources['reference log'] += 1

    known = sorted(cost for cost in estimate.values() if cost is not None)
    if known:
        unknown_cost = known[min(int(UNKNOWN_COST_QUANTILE * len(known)), len(known) - 1)]
    else:
        unknown_cost = float(timeout)

    for input_path, cost in estimate.items():
        estimate[input_path] = (unknown_cost if cost is None else cost) + STARTUP_COST
    return estimate, sources

'''
    group the input scripts into units that have to be run by the same worker

    Input scripts in the same folder may write files that another input script in that
    folder reads, or may overwrite each other's output files.  Those folders are listed
    in the "keep_together" section of the test configuration as patterns relative to the
    top-level examples folder, and all their input scripts become a single unit that is
    run by one worker in the order in which the input scripts are listed.

    return a list of lists of input script paths
'''
def make_work_units(input_list, example_toplevel, keep_together):
    units = []
    grouped = {}
    for input_path in input_list:
        folder = os.path.dirname(example_key(input_path, example_toplevel))
        if any(fnmatch.fnmatch(folder, pattern) for pattern in keep_together):
            if folder not in grouped:
                grouped[folder] = []
                units.append(grouped[folder])
            grouped[folder].append(input_path)
        else:
            units.append([input_path])
    return units

'''
    distribute a list of work units over N workers, most expensive unit first

    This is the "longest processing time first" heuristic, which is also what CTest
    uses to schedule tests with a COST property: the units are sorted by decreasing
    cost and each unit is handed to the worker that has the least work assigned so far.
    The resulting makespan is at worst 4/3 of the optimum minus 1/(3N).

    Each worker's list is ordered by decreasing cost as well, so that a worker that
    gets behind is behind with a short run and not with the longest one.

    return a tuple of (list of N lists of input script paths, list of N estimated
    total times in seconds)
'''
def distribute_by_cost(units, cost, N):
    ordered = sorted(units, key=lambda unit: (-sum(cost[i] for i in unit), unit[0]))

    # heap of (accumulated cost, worker id) with the least loaded worker on top
    workers = [(0.0, idx) for idx in range(N)]
    heapq.heapify(workers)
    sublists = [[] for idx in range(N)]
    loads = [0.0] * N
    for unit in ordered:
        load, idx = heapq.heappop(workers)
        sublists[idx].extend(unit)
        load += sum(cost[i] for i in unit)
        loads[idx] = load
        heapq.heappush(workers, (load, idx))
    return sublists, loads

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
    parser.add_argument("--cost-file", dest="cost_file", default="",
                        help="JSON file with measured run times from a previous run "
                             "(written by merge_results.py --cost-file) used to distribute "
                             "the input scripts evenly over the workers")
    parser.add_argument("--walltime-ref", dest="walltime_ref", default="",
                        help="Reference walltime in seconds used to normalize the reported run times; "
                             "if not given, it is measured by running bench/in.lj")
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
    cost_file = args.cost_file

    # logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode="w")

    # read in the configuration of the tests; it is needed for splitting the input
    # scripts over the workers as well, so it is read before anything else is done
    try:
        with open(configFileName, 'r') as f:
            config = yaml.load(f, Loader=Loader)
    except OSError as err:
        print(f"Cannot read the test configuration file: {err}")
        print("Use --config-file to select one of the config*.yaml files in tools/regression-tests")
        quit()
    print(f"\nRegression test configuration file:\n  {os.path.abspath(configFileName)}")

    # number of procs used for the runs, needed to scale the cost estimates taken
    # from the reference log files; an empty value means "as many as the reference
    # log file with the most procs", which is 4 for most of the examples
    # the reference log files are named after the number of MPI processes times the
    # number of OpenMP threads, so the cost estimates are scaled with that as well
    nprocs_config = config_ntasks(4 if config['nprocs'] == "" else int(config['nprocs']),
                                  config['args'])

    # the runs are killed after this many seconds (must match the default in execute())
    timeout_config = int(config['timeout']) if config.get('timeout', "") != "" else 60

    # folders whose input scripts must be run by the same worker in the listed order
    keep_together = config['keep_together'] if 'keep_together' in config else []

    '''
        split a list of input scripts over the workers and write one input-list-{idx}.txt
        file per worker, ordered so that the most expensive runs are started first
    '''
    def write_input_lists(input_list):
        costs = load_cost_file(cost_file)
        if cost_file:
            print(f"\nRead {len(costs)} measured run times from {cost_file}")
        cost, sources = estimate_costs(input_list, example_toplevel, costs, nprocs_config,
                                       timeout_config)
        units = make_work_units(input_list, example_toplevel, keep_together)
        sublists, loads = distribute_by_cost(units, cost, num_workers)

        total = sum(cost.values())
        msg = (f"\nDistributing {len(input_list)} input scripts in {len(units)} units over "
               f"{num_workers} workers:\n"
               f"  Cost estimates: " + ", ".join(f"{num} from {src}" for src, num in sources.items()) + "\n"
               f"  Estimated total time    : {total:.0f} s\n"
               f"  Estimated time per worker: {min(loads):.0f} s (min), "
               f"{total/num_workers:.0f} s (average), {max(loads):.0f} s (max)")
        print(msg)
        logger.info(msg)

        for idx, sublist in enumerate(sublists):
            with open(f"input-list-{idx}.txt", "w") as f:
                for inp in sublist:
                    f.write(inp + '\n')
        return cost

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

            # distribute the input scripts over the workers by their estimated cost
            write_input_lists(input_list)
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
            # sort the list, since find returns the input scripts in file system order
            input_list = sorted(inp for inp in p.stdout.split('\n') if inp)
            msg = f"\nThere are {len(input_list)} input scripts in total under the {example_toplevel} folder."
            print(msg)
            logger.info(msg)

            # get the list of folders that contain the input scripts
            folder_list = []
            inputs_per_folder = {}
            for input in input_list:
                folder = input.rsplit('/', 1)[0]
                # unique folders in the list
                if folder not in inputs_per_folder:
                    inputs_per_folder[folder] = []
                    folder_list.append(folder)
                inputs_per_folder[folder].append(input)

            # distribute the input scripts over the workers by their estimated cost
            cost = write_input_lists(input_list)

            # distribute the folders over the workers by the summed cost of their inputs
            folder_cost = { folder: sum(cost[inp] for inp in inputs)
                            for folder, inputs in inputs_per_folder.items() }
            sublists, loads = distribute_by_cost([[folder] for folder in folder_list],
                                                 folder_cost, num_workers)
            for idx, sublist in enumerate(sublists):
                with open(f"folder-list-{idx}.txt", "w") as f:
                    for folder in sublist:
                        # each line lists a folder and its number of input scripts
                        f.write(f"{folder} {len(inputs_per_folder[folder])}\n")

            # working on all the folders for now
            example_subfolders = folder_list

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

            # the top-level examples folder is needed to locate bench/in.lj for the
            # reference walltime; derive it from the listed input scripts if it was
            # not given on the command line
            if example_toplevel == "" and len(example_inputs) > 0:
                marker = '/examples/'
                path = os.path.abspath(example_inputs[0])
                if marker in path:
                    example_toplevel = path.split(marker, 1)[0] + '/examples'

        else:
            inplace_input = False

    # if analyze the example folders (and split into separate lists for top-level examples), not running any test
    if analyze == True:
        # measure the reference walltime once and store it, so that the workers can be
        # started with --walltime-ref instead of all benchmarking the LAMMPS binary
        # against each other at the same time
        if os.path.isfile(lmp_binary) and example_toplevel != "":
            walltime_ref = get_reference_walltime(lmp_binary, config, example_toplevel)
            with open("walltime-ref.txt", "w") as f:
                f.write(f"{walltime_ref}\n")
            print(f"\nReference walltime, sec = {walltime_ref} (written to walltime-ref.txt)")
        quit()

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

    # get a reference walltime; when many workers run concurrently, it should be
    # measured once and passed in with --walltime-ref instead, since the workers
    # would otherwise all benchmark against each other at the same time
    if args.walltime_ref:
        walltime_ref = float(args.walltime_ref)
        print(f"\nReference walltime, sec = {walltime_ref}")
    else:
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

        # build the list of work units, each of them a folder and the input scripts
        # to be tested in it.  If a list of input scripts was given, its order is
        # preserved, so that the most expensive runs are started first when the list
        # was written by a run with --analyze; consecutive input scripts from the same
        # folder are collected into one work unit.  Otherwise all input scripts in a
        # folder are tested in alphabetical order.
        work_units = []
        if len(example_inputs) > 0:
            for input in example_inputs:
                folder, name = input.rsplit('/', 1)
                if work_units and work_units[-1][0] == folder:
                    work_units[-1][1].append(name)
                else:
                    work_units.append((folder, [name]))
        else:
            for directory in example_subfolders:
                work_units.append((directory, None))

        for directory, input_list in work_units:

            if os.path.exists(directory) is False:
                continue

            # change to the directory where the input script and data files are located
            print("-"*80)
            print("Entering " + directory)
            logger.info("Entering " + directory)
            os.chdir(directory)

            # test all input scripts in the folder if none were selected explicitly
            if input_list is None:
                input_list = sorted(glob.glob("in.*"))
            else:
                input_list = [inp for inp in input_list if os.path.isfile(inp)]

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
