'''
UPDATE: Feb 8, 2024:
  Launching the LAMMPS binary under testing using a configuration defined in a yaml file (e.g. config.yaml).
  Comparing the output thermo with that in the existing log file (with the same nprocs)
    + data in the log files are extracted and converted into yaml data structure
    + using the in place input scripts, no need to add REG markers to the input scripts
  This way we can:
    + launch tests with mpirun with multiple procs
    + specify what LAMMPS binary version to test (e.g., testing separate builds)
    + simplify the build configuration (no need to build the Python module)
  NOTE: Need to allow to tolerances specified for invidual input scripts,
        or each config.yaml is for a set of example folders

  Example usage:
    1) Simple use (using the provided tools/regression-tests/config.yaml and the examples/ folder at the top level)
       python3 run_tests.py --lmp-bin=/path/to/lmp_binary
    2) Use a custom testing configuration
       python3 run_tests.py --lmp-bin=/path/to/lmp_binary --config-file=/path/to/config/file/config.yaml
    3) Specify a list of example folders with a modifed configuration (e.g. different tolerances)
       python3 run_tests.py --lmp-bin=/path/to/lmp_binary \
          --example-folders="/path/to/examples/folder1;/path/to/examples/folder2" \
          --config-file=/path/to/config/file/config.yaml 
'''

import os
import datetime
import fnmatch
import subprocess
from argparse import ArgumentParser

# need "pip install pyyaml numpy"
import yaml
import numpy as np

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
  inputFileName:  input file with comments #REG:ADD and #REG:SUB as markers
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
  yamlFileName: input YAML file with thermo structured
    as described in https://docs.lammps.org/Howto_structured_data.html
  return: thermo, which is a list containing a dictionary for each run
    where the tag "keywords" maps to the list of thermo header strings
    and the tag “data” has a list of lists where the outer list represents the lines
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
  Convert an existing log.lammps file into a thermo yaml style log
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
    for line in lines:
      if "Step" in line:
        line.strip()
        keywords = line.split()
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
  return the list of installed packages
'''
def get_lammps_build_configuration(lmp_binary):
  cmd_str = lmp_binary + " -h"
  p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
  output = p.stdout.split('\n')
  packages = ""
  reading = False
  row = 0
  for l in output:
    if l != "":
      if l == "Installed packages:":
        reading = True
        n = row
      if "List of individual style options" in l:
        reading = False
      if reading == True and row > n:
        packages += l.strip() + " "

    if "OS:" in l:
      operating_system = l
    if "Git info" in l:
      GitInfo = l
 
    row += 1

  packages = packages.strip()

  row = 0
  compile_flags = ""
  for l in output:
    if l != "":
      if "-DLAMMPS" in l:
        compile_flags += " " + l.strip()
 
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

  #output = p.stdout.split('\n')
  output = p.stdout
  # process output to handle failed runs

  return cmd_str, output


'''
  attempt to plug in the REG markers before each run command
  #REG:ADD thermo 10
  #REG:ADD thermo_style yaml
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

def has_markers(input):
  with open(input) as f:
    if '#REG' in f.read():
      return True
  return False

'''
  Iterate over a list of input files using the testing configuration
  return total number of tests, and the number of tests with failures
'''
def iterate(input_list, config, results, removeAnnotatedInput=False):
  EPSILON = np.float64(config['epsilon'])
  nugget = float(config['nugget'])

  num_tests = len(input_list)
  num_passed = 0
  test_id = 0

  # using REG-commented input scripts
  using_markers = False

  # iterate over the input scripts
  for input in input_list:

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

    # check if a log file exists in the current folder: log.[date].basename.[nprocs]
    basename = input_test.replace('in.','')    
    logfile_exist = False

    # if there are multiple log files for different number of procs, pick the maximum number
    pattern = f'log.*.{basename}.*'
    max_np = 1
    for file in os.listdir('.'):
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
      print(f"Cannot find reference log file with {pattern}.")
      # try to read in the thermo yaml output from the working directory
      thermo_ref_file = 'thermo.' + input + '.yaml'
      file_exist = os.path.isfile(thermo_ref_file)
      if file_exist == True:
        thermo_ref = extract_thermo(thermo_ref_file)
        num_runs_ref = len(thermo_ref)
      else:
        print(f"SKIPPED: {thermo_ref_file} does not exist in the working directory.")
        result.status = "skipped"
        results.append(result)
        continue

    # using the LAMMPS python module (for single-proc runs)
    #  lmp = lammps()
    #  lmp.file(input_test)
    
    # or more customizable with config.yaml
    cmd_str, output = execute(lmp_binary, config, input_test)

    # restore the nprocs value in the configuration
    config['nprocs'] = saved_nprocs

    # process thermo output
    thermo = extract_data_to_yaml("log.lammps")

    num_runs = len(thermo)
    if num_runs == 0:
      print(f"ERROR: Failed with the running with {input_test}. The run terminated with the following output:\n")
      print(f"{output}")
      result.status = "error"
      results.append(result)
      continue 

    print(f"Comparing thermo output from log.lammps with the reference log file {thermo_ref_file}")
    if num_runs != num_runs_ref:
      print(f"ERROR: Number of runs in log.lammps ({num_runs}) is not the same as that in the reference log ({num_runs_ref})")
      result.status = "error"
      results.append(result)
      continue

    # comparing output vs reference values
    width = 20
    if verbose == True:
      print("Quantities".ljust(width) + "Output".center(width) + "Reference".center(width) + "Abs Diff Check".center(width) +  "Rel Diff Check".center(width))
    
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
          print(f"{thermo[irun]['keywords'][i].ljust(width)} {str(val).rjust(20)} {str(ref).rjust(20)} {abs_diff_check.rjust(20)} {rel_diff_check.rjust(20)}")

    if num_abs_failed > 0:
      print(f"{num_abs_failed} absolute diff checks failed with the specified tolerances.")
      result.status = "failed"
      for i in failed_abs_output:
        print(f"- {i}")
    if num_rel_failed > 0:
      print(f"{num_rel_failed} relative diff checks failed with the specified tolerances.")
      result.status = "failed"
      for i in failed_rel_output:
        print(f"- {i}")
    if num_abs_failed == 0 and num_rel_failed == 0:
      print(f"All {num_checks} checks passed.")
      result.status = "passed"        
      num_passed = num_passed + 1

    results.append(result)

    print("-"*(5*width+4))
    test_id = test_id + 1

    # remove the annotated input script
    if removeAnnotatedInput == True:
      cmd_str = "rm " + input_test
      os.system(cmd_str)

  return num_passed

'''
  TODO:
    - automate annotating the example input scripts if thermo style is multi (e.g. examples/peptide)
'''
if __name__ == "__main__":

  # default values
  lmp_binary = ""
  configFileName = "config.yaml"
  example_subfolders = []
  genref = False
  verbose = False
  output_file = "output.xml"

  # parse the arguments
  parser = ArgumentParser()
  parser.add_argument("--lmp-bin", dest="lmp_binary", default="", help="LAMMPS binary")
  parser.add_argument("--config-file", dest="config_file", default="config.yaml",
                      help="Configuration YAML file")
  parser.add_argument("--example-folders", dest="example_folders", default="", help="Example subfolders")
  parser.add_argument("--gen-ref",dest="genref", action='store_true', default=False,
                      help="Generating reference data")
  parser.add_argument("--verbose",dest="verbose", action='store_true', default=False,
                      help="Verbose output")
  parser.add_argument("--output",dest="output", default="output.xml", help="Output file")

  args = parser.parse_args()

  lmp_binary = os.path.abspath(args.lmp_binary)
  configFileName = args.config_file
  output_file = args.output
  if args.example_folders != "":
    example_subfolders = args.example_folders.split(';')
    print("Example folders:")
    print(example_subfolders)
  genref = args.genref
  verbose = args.verbose

  # read in the configuration of the tests
  with open(configFileName, 'r') as f:
    config = yaml.load(f, Loader=Loader)
    absolute_path = os.path.abspath(configFileName)
    print(f"Regression tests with settings defined in {absolute_path}")
 
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
  print(f"- List of installed packages: {packages}")
  
  # Using in place input scripts
  inplace_input = True
  test_cases = []

  # if the example folders are not specified from the command-line argument -example-folders
  if len(example_subfolders) == 0:
    example_subfolders.append("../../examples/melt")
    example_subfolders.append('../../examples/flow')
    example_subfolders.append('../../examples/indent')
    example_subfolders.append('../../examples/shear')
    example_subfolders.append('../../examples/steinhardt')

    # prd  log file parsing issue
    # neb  log file parsing issue
    # snap log files obsolete?

    # append the example subfolders depending on the installed packages
    if 'ASPHERE' in packages:
      #example_subfolders.append('../../examples/ASPHERE/ellipsoid')
      example_subfolders.append('../../examples/ellipse')

    if 'CORESHELL' in packages:
      example_subfolders.append('../../examples/coreshell')

    if 'MOLECULE' in packages:
      example_subfolders.append('../../examples/micelle')
      # peptide thermo_style as multi
      #example_subfolders.append('../../examples/peptide')

    if 'GRANULAR' in packages:
      example_subfolders.append('../../examples/granular')
      example_subfolders.append('../../examples/pour')

    if 'AMOEBA' in packages:
      example_subfolders.append('../../examples/amoeba')

    if 'BODY' in packages:
      example_subfolders.append('../../examples/body')

    if 'BPM' in packages:
      example_subfolders.append('../../examples/bpm/impact')
      example_subfolders.append('../../examples/bpm/pour')

    if 'COLLOID' in packages:
      example_subfolders.append('../../examples/colloid')

    if 'CRACK' in packages:
      example_subfolders.append('../../examples/crack')

    if 'DIELECTRIC' in packages:
      example_subfolders.append('../../examples/PACKAGES/dielectric')

    if 'DIPOLE' in packages:
      example_subfolders.append('../../examples/dipole')


    if 'DPD-BASIC' in packages:
      example_subfolders.append('../../examples/PACKAGES/dpd-basic/dpd')
      example_subfolders.append('../../examples/PACKAGES/dpd-basic/dpdext')
      example_subfolders.append('../../examples/PACKAGES/dpd-basic/dpd_tstat')
      example_subfolders.append('../../examples/PACKAGES/dpd-basic/dpdext_tstat')

    if 'MANYBODY' in packages:
      example_subfolders.append('../../examples/tersoff')
      example_subfolders.append('../../examples/vashishta')
      example_subfolders.append('../../examples/threebody')

    if 'RIGID' in packages:
      example_subfolders.append('../../examples/rigid')

    if 'SRD' in packages:
      example_subfolders.append('../../examples/srd')
    
  all_results = []
  if inplace_input == True:

    # save current working dir
    p = subprocess.run("pwd", shell=True, text=True, capture_output=True)
    pwd = p.stdout.split('\n')[0]
    pwd = os.path.abspath(pwd)
    print("Working directory: " + pwd)

    # change dir to a folder under examples/, need to use os.chdir()
    # TODO: loop through the subfolders under examples/, depending on the installed packages

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

      print(f"List of input scripts: {input_list}")
      total_tests += len(input_list)

      # iterate through the input scripts
      results = []
      num_passed = iterate(input_list, config, results)
      passed_tests += num_passed

      all_results.extend(results)

      # get back to the working dir
      os.chdir(pwd)

  else:
    # or using the input scripts in the working directory -- for debugging purposes
    input_list=['in.lj', 'in.rhodo', 'in.eam']
    total_tests = len(input_list)
    results = []
    passed_tests = iterate(input_list, config, results)

  print("Summary:")
  print(f" - {passed_tests} passed / {total_tests} tests")
  print(f" - Details are given in {output_file}.")

  # generate a JUnit XML file
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