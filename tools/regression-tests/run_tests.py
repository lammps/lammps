#!/usr/bin/env python3
'''
UPDATE: July 21, 2024:
  Launching the LAMMPS binary under testing using a configuration defined in a yaml file (e.g. config.yaml).
  Comparing the output thermo with that in the existing log file (with the same nprocs)
    + data in the log files are extracted and converted into yaml data structure
    + using the in place input scripts, no need to add REG markers to the input scripts

With the current features, users can:
    + specify which LAMMPS binary version to test (e.g., the version from a commit, or those from `lammps-testing`)
    + specify the examples subfolders (thus the reference log files) seperately (e.g. from other LAMMPS versions or commits)
    + specify tolerances for individual quantities for any input script to override the global values
    + launch tests with `mpirun` with all supported command line features (multiple procs, multiple paritions, and suffices)
    + skip certain input files if not interested, or no reference log file exists
    + simplify the main LAMMPS builds, as long as a LAMMPS binary is available

Limitations:
    - input scripts use thermo style multi (e.g., examples/peptide) do not work with the expected thermo output format
    - input scripts that require partition runs (e.g. examples/neb) need a separate config file, e.g. "args: --partition 2x1"
    - testing accelerator packages (GPU, INTEL, KOKKOS, OPENMP) need separate config files, "args: -sf omp -pk omp 4"

TODO:
    + keep track of the testing progress to resume the testing from the last checkpoint
    + distribute the input list across multiple processes via multiprocessing, or
      split the list of input scripts into separate runs (there are 800+ input script under the top-level examples)
    + be able to be invoked from run_tests in the lammps-testing infrastruture

The following Python packages need to be installed into an activated environment:
    
    python3 -m venv testing-env
    source testing-env/bin/activate
    pip install numpy pyyaml junit_xml

Example usage:
    1) Simple use (using the provided tools/regression-tests/config.yaml and the examples/ folder at the top level)
       python3 run_tests.py --lmp-bin=/path/to/lmp_binary
    2) Use a custom testing configuration
       python3 run_tests.py --lmp-bin=/path/to/lmp_binary --config-file=/path/to/config/file/config.yaml
    3) Specify a list of example folders with a modifed configuration (e.g. different tolerances)
       python3 run_tests.py --lmp-bin=/path/to/lmp_binary \
          --example-folders="/path/to/examples/folder1;/path/to/examples/folder2" \
          --config-file=/path/to/config/file/config.yaml 
    4) Test a LAMMPS binary with the whole top-level /examples folder in a LAMMPS source tree
       python3 run_tests.py --lmp-bin=/path/to/lmp_binary --example-top-level=/path/to/lammps/examples
'''

import os
import sys
import datetime
import fnmatch
import subprocess
from argparse import ArgumentParser

from multiprocessing import Pool

import logging
# need "pip install numpy pyyaml"
import numpy as np
import yaml

# need "pip install junit_xml"
from junit_xml import TestSuite, TestCase

try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader

class TestResult:
  def __init__(self, name, output=None, time=None, checks=0, status=None):
    self.name = name
    self.output = output
    self.time = time
    self.checks = 0
    self.status = status


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
            if "Step" in line:
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
    try:
        yaml_struct = yaml.load_all(docs, Loader=Loader)
        thermo = list(yaml_struct)
    except yaml.YAMLError as exc:
        if hasattr(exc, 'problem_mark'):
            mark = exc.problem_mark
            print(f"Error parsing {inputFileName} at line {mark.line}, column {mark.column+1}.")
        else:
            print (f"Something went wrong while parsing {inputFileName}.")
            print(docs)
    return thermo

'''
    return a tuple of the list of installed packages, OS, GitInfo and compile_flags
'''
def get_lammps_build_configuration(lmp_binary):
    cmd_str = lmp_binary + " -h"
    p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
    output = p.stdout.split('\n')
    packages = ""
    reading = False
    operating_system = ""
    GitInfo = ""
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
  
        row += 1

    packages = packages.strip()

    row = 0
    compile_flags = ""
    for line in output:
        if line != "":
            if "-DLAMMPS" in line:
                compile_flags += " " + line.strip()
  
        row += 1

    return packages.split(" "), operating_system, GitInfo, compile_flags

'''
    launch LAMMPS using the configuration defined in the dictionary config with an input file
    TODO:
        - generate new reference values if needed
        - wrap subprocess with try/catch to handle exceptions
'''
def execute(lmp_binary, config, input_file_name, generate_ref_yaml=False):
    cmd_str = config['mpiexec'] + " " + config['mpiexec_numproc_flag'] + " " + config['nprocs'] + " "
    cmd_str += lmp_binary + " -in " + input_file_name + " " + config['args']
    print(f"Executing: {cmd_str}")
    p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)

    return cmd_str, p.stdout, p.stderr, p.returncode

'''
    split a list into a list of N sublists
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
    Iterate over a list of input files using the given lmp_binary, the testing configuration
    return test results, as a list of TestResult instances

    To map a function to individual workers:

    def func(input1, input2, output):
        # do smth
        return result

    # args is a list of num_workers tuples, each tuple contains the arguments passed to the function executed by a worker
    args = []
    for i in range(num_workers):
        args.append((input1, input2, output))

    with Pool(num_workers) as pool:   
        results = pool.starmap(func, args)

'''
def iterate(lmp_binary, input_list, config, results, removeAnnotatedInput=False, output=None):
    EPSILON = np.float64(config['epsilon'])
    nugget = float(config['nugget'])

    num_tests = len(input_list)
    num_passed = 0
    test_id = 0

    # using REG-commented input scripts, now turned off (False)
    using_markers = False

    # iterate over the input scripts
    for input in input_list:

        # skip the input file if listed
        if 'skip' in config:
            if input in config['skip']:
                logger.info(f"SKIPPED: {input} as specified in the configuration file {configFileName}")
                test_id = test_id + 1
                continue

        str_t = "\nRunning " + input + f" ({test_id+1}/{num_tests})"

        result = TestResult(name=input, output="", time="", status="passed")

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
                
                str_t = "\nRunning " + input_test + f" ({test_id+1}/{num_tests})"
        else:
            input_test = input

        print(str_t)
        print(f"-"*len(str_t))
        logger.info(str_t)
        logger.info(f"-"*len(str_t))

        # check if a log file exists in the current folder: log.DDMMMYY.basename.[nprocs]
        basename = input_test.replace('in.','')
        logfile_exist = False

        # if there are multiple log files for different number of procs, pick the maximum number
        cmd_str = "ls log.*"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        logfile_list = p.stdout.split('\n')
        logfile_list.remove('')

        max_np = 1
        for file in logfile_list:
            # looks for pattern log.{date}.{basename}.g++.{nprocs}
            # get the date from the log files
            date = file.split('.',2)[1]
            pattern = f'log.{date}.{basename}.*'
            if fnmatch.fnmatch(file, pattern):
                p = file.rsplit('.', 1)
                if max_np < int(p[1]):
                    max_np = int(p[1])
                    logfile_exist = True
                    thermo_ref_file = file

        # if the maximum number of procs is different from the value in the configuration file
        #      then override the setting for this input script
        saved_nprocs = config['nprocs']
        if max_np != int(config['nprocs']):
            config['nprocs'] = str(max_np)

        if logfile_exist:
            thermo_ref = extract_data_to_yaml(thermo_ref_file)
            num_runs_ref = len(thermo_ref)
        else:
            logger.info(f"Cannot find a reference log file {thermo_ref_file} for {input_test}.")
            # try to read in the thermo yaml output from the working directory
            thermo_ref_file = 'thermo.' + input + '.yaml'
            file_exist = os.path.isfile(thermo_ref_file)
            if file_exist == True:
                thermo_ref = extract_thermo(thermo_ref_file)
                num_runs_ref = len(thermo_ref)
            else:
                logger.info(f"SKIPPED: {thermo_ref_file} does not exist in the working directory.")
                result.status = "skipped due to missing the log file"
                results.append(result)
                test_id = test_id + 1
                continue

        # or more customizable with config.yaml
        cmd_str, output, error, returncode = execute(lmp_binary, config, input_test)

        # restore the nprocs value in the configuration
        config['nprocs'] = saved_nprocs

        # process error code from the run
        if os.path.isfile("log.lammps") == False:
            logger.info(f"ERROR: No log.lammps generated with {input_test} with return code {returncode}. Check the {log_file} for the run output.\n")
            logger.info(f"\n{input_test}:")
            logger.info(f"\n{error}")
            test_id = test_id + 1
            continue

        # process thermo output from the run
        thermo = extract_data_to_yaml("log.lammps")

        num_runs = len(thermo)
        if num_runs == 0:
            logger.info(f"The run terminated with {input_test} gives the following output:\n")
            logger.info(f"\n{output}")
            if "Unrecognized" in output:
                result.status = "unrecognized command"
            else:
                result.status = "error"
                logger.info(f"ERROR: Failed with {input_test} due to {result.status}.\n")
                results.append(result)
                test_id = test_id + 1
                continue

        print(f"Comparing thermo output from log.lammps against the reference log file {thermo_ref_file}")
        if num_runs != num_runs_ref:
            logger.info(f"ERROR: Number of runs in log.lammps ({num_runs}) is not the same as that in the reference log ({num_runs_ref})")
            result.status = "error"
            results.append(result)
            test_id = test_id + 1
            continue

        # comparing output vs reference values
        width = 20
        if verbose == True:
            print("Quantities".ljust(width) + "Output".center(width) + "Reference".center(width) +
                "Abs Diff Check".center(width) +  "Rel Diff Check".center(width))
        
        # check if overrides for this input scipt is specified
        overrides = {}
        if 'overrides' in config:
            if input_test in config['overrides']:
                overrides = config['overrides'][input_test]

        # iterate through all num_runs

        num_abs_failed = 0
        num_rel_failed = 0
        failed_abs_output = []
        failed_rel_output = []
        num_checks = 0

        for irun in range(num_runs):
            num_fields = len(thermo[irun]['keywords'])

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
                print(f"{thermo[irun]['keywords'][i].ljust(width)} {str(val).rjust(20)} {str(ref).rjust(20)} "
                    "{abs_diff_check.rjust(20)} {rel_diff_check.rjust(20)}")

        if num_abs_failed > 0:
            print(f"{num_abs_failed} absolute diff checks failed with the specified tolerances.")
            result.status = "failed"
            if verbose == True:
                for i in failed_abs_output:
                    print(f"- {i}")
        if num_rel_failed > 0:
            print(f"{num_rel_failed} relative diff checks failed with the specified tolerances.")
            result.status = "failed"
            if verbose == True:
                for i in failed_rel_output:
                    print(f"- {i}")
        if num_abs_failed == 0 and num_rel_failed == 0:
            print(f"All {num_checks} checks passed.")
            result.status = "passed"        
            num_passed = num_passed + 1

        results.append(result)

        str_t = f"Finished " + input_test
        print(str_t)
        print("-"*(5*width+4))
        test_id = test_id + 1

        # remove the annotated input script
        if removeAnnotatedInput == True:
            cmd_str = "rm " + input_test
            os.system(cmd_str)

    return num_passed

'''
  TODO:


'''
if __name__ == "__main__":

    # default values
    lmp_binary = ""
    configFileName = "config.yaml"
    example_subfolders = []
    example_toplevel = ""
    genref = False
    verbose = False
    output_file = "output.xml"
    log_file = "run.log"
    dry_run = False
    
    # distribute the total number of input scripts over the workers
    num_workers = 1

    # parse the arguments
    parser = ArgumentParser()
    parser.add_argument("--lmp-bin", dest="lmp_binary", default="", help="LAMMPS binary")
    parser.add_argument("--config-file", dest="config_file", default=configFileName,
                        help="Configuration YAML file")
    parser.add_argument("--example-top-level", dest="example_toplevel", default="", help="Example top-level")
    parser.add_argument("--example-folders", dest="example_folders", default="", help="Example subfolders")
    parser.add_argument("--num-workers", dest="num_workers", default=1, help="Number of workers")
    parser.add_argument("--gen-ref",dest="genref", action='store_true', default=False,
                        help="Generating reference data")
    parser.add_argument("--verbose",dest="verbose", action='store_true', default=False,
                        help="Verbose output")
    parser.add_argument("--output",dest="output", default=output_file, help="Output file")
    parser.add_argument("--logfile",dest="logfile", default=log_file, help="Log file")
    parser.add_argument("--dry-run",dest="dry_run", action='store_true', default=False,
                        help="Only report statistics")

    args = parser.parse_args()

    lmp_binary = os.path.abspath(args.lmp_binary)
    configFileName = args.config_file
    output_file = args.output
    num_workers = args.num_workers

    # example_toplevel is where all the examples subfolders reside
    if args.example_toplevel != "":
       example_toplevel = args.example_toplevel
    if args.example_folders != "":
        example_subfolders = args.example_folders.split(';')
       
    genref = args.genref
    verbose = args.verbose
    log_file = args.logfile
    dry_run = args.dry_run

    # logging
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=log_file, level=logging.INFO, filemode="w")

    # read in the configuration of the tests
    with open(configFileName, 'r') as f:
        config = yaml.load(f, Loader=Loader)
        absolute_path = os.path.abspath(configFileName)
        print(f"Regression tests with settings defined in the configuration file {absolute_path}")
  
    # check if lmp_binary is specified in the config yaml
    if lmp_binary == "":
        if config['lmp_binary'] == "":
            print("Needs a valid LAMMPS binary")
            quit()
        else:
            lmp_binary = os.path.abspath(config['lmp_binary'])

    # print out the binary info
    packages, operating_system, GitInfo, compile_flags = get_lammps_build_configuration(lmp_binary)
    print("LAMMPS build info:")
    print(f"- {operating_system}")
    print(f"- {GitInfo}")
    print(f"- Active compile flags: {compile_flags}")
    print(f"- List of {len(packages)} installed packages: {packages}")

    if len(example_subfolders) > 0:
        print("- Example folders to test:")
        print(example_subfolders)
    if example_toplevel != "":
        print("- Top-level example folder:")
        print(example_toplevel)

    folder_list = []
    if len(example_toplevel) != 0:
        # getting the list of all the input files because there are subfolders (e.g. PACKAGES) under the top level
        cmd_str = f"find {example_toplevel} -name \"in.*\" "
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        input_list = p.stdout.split('\n')
        input_list.remove("")

        # find out which folder to cd into to run the input script
        for input in input_list:
            folder = input.rsplit('/', 1)[0]
            folder_list.append(folder)
        print(f"There are {len(input_list)} input scripts in total under the {example_toplevel} folder.")

        # divide the list of input scripts into num_workers chunks
        sublists = divide_into_N(input_list, num_workers)

    # if only statistics, not running anything
    if dry_run == True:
        quit()

    # Using in place input scripts
    inplace_input = True
    test_cases = []

    # if the example folders are not specified from the command-line argument --example-folders
    # then use the path from --example-top-folder
    if len(example_subfolders) == 0:

        # get the input file list, for now the first in the sublist
        # TODO: generate a list of tuples, each tuple contains a folder list for a worker,
        #       then use multiprocessing.Pool starmap()
        folder_list = []
        for input in sublists[0]:
            folder = input.rsplit('/', 1)[0]
            # unique folders in the list
            if folder not in folder_list:
                folder_list.append(folder)

        example_subfolders = folder_list

    all_results = []

    # default setting is to use inplace_input
    if inplace_input == True:
        # save current working dir
        p = subprocess.run("pwd", shell=True, text=True, capture_output=True)
        pwd = p.stdout.split('\n')[0]
        pwd = os.path.abspath(pwd)
        print("Working directory: " + pwd)

        # change dir to a folder under examples/, need to use os.chdir()
        # TODO: loop through the subfolders under examples/, depending on the installed packages

        '''
        args = []
        for i in range(num_workers):
            args.append((input1, input2, output))

        with Pool(num_workers) as pool:   
            results = pool.starmap(func, args)
        '''
        total_tests = 0
        passed_tests = 0

        for directory in example_subfolders:

            p = subprocess.run("pwd", shell=True, text=True, capture_output=True)
            print("\nEntering " + directory)
            os.chdir(directory)

            cmd_str = "ls in.*"
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list.remove('')

            print(f"List of {len(input_list)} input scripts: {input_list}")
            total_tests += len(input_list)

            # iterate through the input scripts
            results = []
            num_passed = iterate(lmp_binary, input_list, config, results)
            passed_tests += num_passed

            # append the results to the all_results list
            all_results.extend(results)

            # get back to the working dir
            os.chdir(pwd)

    else:
        # or using the input scripts in the working directory -- for debugging purposes
        input_list=['in.lj']
        total_tests = len(input_list)
        results = []
        passed_tests = iterate(input_list, config, results)

        all_results.extend(results)

    # print out summary
    print("Summary:")
    print(f" - {passed_tests} numerical tests passed / {total_tests} tests")
    print(f" - Details are given in {output_file}.")

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