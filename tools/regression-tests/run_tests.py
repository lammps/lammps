'''
UPDATE: Jan 30, 2024:
  Launching the LAMMPS binary under testing using a configuration defined in a yaml file (e.g. config.yaml).
  Comparing the output thermo with that in the existing log file (with the same nprocs)
    + data in the log files are extracted and converted into yaml data structure
    + using the in place input scripts, no need to add REG markers to the input scripts
  This way we can:
    + launch tests with mpirun with multiple procs
    + specify what LAMMPS binary version to test 
    + simplify the build configuration
  Example usage:
    1) Simple use:
         python3 run_tests.py --lmp-bin=/path/to/lmp_binary
    2) Use a custom testing configuration
         python3 run_tests.py --lmp-bin=/path/to/lmp_binary --config-file=/path/to/config/file/config.yaml


--------------------------------------------------------------------------------------------------------------
Original plan: using the LAMMPS Python module
  The following steps are for setting up regression tests with the LAMMPS Python module

  0) Create a virtual environment, and activate it

    python -m venv lmp-venv
    source lmp-venv/bin/activate
    PYTHON_EXECUTABLE=`which python`
    INSTALL_PREFIX=$(${PYTHON_EXECUTABLE} -c "import sys; print(sys.prefix)")

  1) Build LAMMPS as a shared lib and install the LAMMPS python module into the virtual environment

    git clone https://github.com/lammps/lammps.git lammps
    cd lammps
    mkdir build && cd build
    cmake ../cmake/presets/basic.cmake -DBUILD_SHARED_LIBS=on -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX
    make -j4
    make install-python

  2) Run this script, after activating the virtual environment and having the input scripts with markers ready:

    source lmp-venv/bin/activate
    python3 run_tests.py

'''

#from lammps import lammps
import os
import fnmatch
import re, yaml
import numpy as np
import subprocess

from argparse import ArgumentParser

try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader

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
  yamlFileName: input YAML file with thermo structured as described in https://docs.lammps.org/Howto_structured_data.html
  return: thermo, which is a list containing a dictionary for each run where the tag "keywords" maps to the list
    of thermo header strings and the tag “data” has a list of lists where the outer list represents the lines
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
        data = line.split()
        docs += " - ["
        for field in enumerate(data):
          docs += field[1] + ", "
        docs += "]\n"

    # load the docs into a YAML data struture
    #print(docs)
    thermo = list(yaml.load_all(docs, Loader=Loader))
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
  output = p.stdout.split('\n')


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
'''
def iterate(input_list, config, removeAnnotatedInput=False):
  EPSILON = np.float64(config['epsilon'])
  nugget = float(config['nugget'])

  nprocs = int(config['nprocs'])

  num_tests = len(input_list)
  test_id = 0
  using_markers = False
  
  # iterate over the input scripts
  for input in input_list:

    str_t = "\nRunning " + input + f" ({test_id+1}/{num_tests})"

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

    # using the LAMMPS python module (for single-proc runs)
    #  lmp = lammps()
    #  lmp.file(input_test)
    
    # or more customizable with config.yaml
    execute(lmp_binary, config, input_test)

    # process thermo output
    #thermo = extract_thermo("log.lammps")
    thermo = extract_data_to_yaml("log.lammps")

    num_runs = len(thermo)
    if num_runs == 0:
      print(f"Failed with {input_test}. Check if the run with this input script completed normally.\n")
      continue 

    # check if a log file exists in the current folder: log.[date].basename.[nprocs]
    basename = input_test.replace('in.','')    
    logfile_exist = False
    pattern = f'log.*.{basename}.*.{nprocs}'
    for file in os.listdir('.'):
      if fnmatch.fnmatch(file, pattern):
        logfile_exist = True
        break

    if logfile_exist:
      thermo_ref = extract_data_to_yaml(file)
      num_runs_ref = len(thermo_ref)
    else:
      print(f"Cannot find reference log file with {pattern}")
      # read in the thermo yaml output from the working directory
      thermo_ref_file = 'thermo.' + input + '.yaml'
      file_exist = os.path.isfile(thermo_ref_file)
      if file_exist == True:
        thermo_ref = extract_thermo(thermo_ref_file)
      else:
        print(f"SKIPPED: {thermo_ref_file} does not exist")
        continue

    print(f"Comparing thermo output from log.lammps with the reference log {file}")
    if num_runs != num_runs_ref:
      print(f"ERROR: Number of runs in log.lammps ({num_runs}) is not the same as that in the reference log ({num_runs_ref})")
      continue
    
    # comparing output vs reference values
    width = 20
    print("Quantities".ljust(width) + "Output".center(width) + "Reference".center(width) + "Abs Diff Check".center(width) +  "Rel Diff Check".center(width))
    irun = 0
    num_fields = len(thermo[irun]['keywords'])
    
    # get the total number of the thermo output lines
    nthermo_steps = len(thermo[irun]['data'])
    # get the output at the last timestep
    thermo_step = nthermo_steps - 1
    #print(f"nthermo_steps = {nthermo_steps}")
    num_abs_failed = 0
    num_rel_failed = 0
    for i in range(num_fields):
      quantity = thermo[0]['keywords'][i]

      val = thermo[irun]['data'][thermo_step][i]
      ref = thermo_ref[irun]['data'][thermo_step][i]
      abs_diff = abs(float(val) - float(ref))
      
      if abs(float(ref)) > EPSILON:
        rel_diff = abs(float(val) - float(ref))/abs(float(ref))
      else:
        rel_diff = abs(float(val) - float(ref))/abs(float(ref)+nugget)
      

      abs_diff_check = "PASSED"
      rel_diff_check = "PASSED"

      if quantity in config['tolerance']:
        if abs_diff > float(config['tolerance'][quantity]['abs']):
          abs_diff_check = "FAILED"
          num_abs_failed = num_abs_failed + 1
        if rel_diff > float(config['tolerance'][quantity]['rel']):
          rel_diff_check = "FAILED"
          num_rel_failed = num_rel_failed + 1
      else:
        abs_diff_check = "N/A"
        rel_diff_check = "N/A"

      print(f"{thermo[irun]['keywords'][i].ljust(width)} {str(val).rjust(20)} {str(ref).rjust(20)} {abs_diff_check.rjust(20)} {rel_diff_check.rjust(20)}")
    if num_abs_failed > 0:
      print(f"{num_abs_failed} abs checks failed")
    elif num_rel_failed > 0:
      print(f"{num_rel_failed} rel checks failed")
    else:
      print("All checks passed. (N/A means tolerance not defined in the config file.)")
    

    print("-"*(5*width+4))
    test_id = test_id + 1

    # remove the annotated input script
    if removeAnnotatedInput == True:
      cmd_str = "rm " + input_test
      os.system(cmd_str)


'''
  TODO:
    - automate annotating the example input scripts of the installed packages
'''
if __name__ == "__main__":

  # default values
  lmp_binary = ""
  configFileName = "config.yaml"
  genref = False

  # parse the arguments
  parser = ArgumentParser()
  parser.add_argument("--lmp-bin", dest="lmp_binary", default="", help="LAMMPS binary")
  parser.add_argument("--config-file", dest="config_file", default="config.yaml",
                    help="Configuration YAML file")
  parser.add_argument("--gen-ref",dest="genref", action='store_true', default=False,
                    help="Generating reference data")

  args = parser.parse_args()

  lmp_binary = os.path.abspath(args.lmp_binary)
  config_file= args.config_file
  genref = args.genref
  
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
       lmp_binary = config['lmp_binary']

  # print out the binary info
  packages, operating_system, GitInfo, compile_flags = get_lammps_build_configuration(lmp_binary)
  print("LAMMPS build info:")
  print(f"- {operating_system}")
  print(f"- {GitInfo}")
  print(f"- Active compile flags: {compile_flags}")
  print(f"- List of installed packages: {packages}")
  
  # Using in place input scripts
  inplace_input = True

  example_subfolders = []
  example_subfolders.append("../../examples/melt")

  # append the example subfolders depending on the installed packages
  if 'MOLECULE' in packages:
    example_subfolders.append('../../examples/micelle')

  #if 'ASPHERE' in packages:
  #  example_subfolders.append('../../examples/ASPHERE/ellipsoid')

  if 'COLLOID' in packages:
    example_subfolders.append('../../examples/colloid')

  if 'DIPOLE' in packages:
    example_subfolders.append('../../examples/dipole')

  if 'MANYBODY' in packages:
    example_subfolders.append('../../examples/tersoff')

  if inplace_input == True:

    # save current working dir
    p = subprocess.run("pwd", shell=True, text=True, capture_output=True)
    pwd = p.stdout.split('\n')[0]
    print("Working directory: " + pwd)

    # change dir to a folder under examples/
    # TODO: loop through the subfolders under examples/, depending on the installed packages

    for directory in example_subfolders:

      print("\nEntering " + directory)
      os.chdir(directory)

      # create a symbolic link to the lammps binary at the present directory
      cmd_str = "ln -s " + lmp_binary + " lmp"
      os.system(cmd_str)

      cmd_str = "ls in.*"
      p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
      input_list = p.stdout.split('\n')
      input_list.remove('')

      print("List of input scripts:")
      print(input_list)

      # iterate through the input scripts
      iterate(input_list, config)

      # unlink the symbolic link
      cmd_str = "unlink lmp"
      os.system(cmd_str)
      # get back to the working dir
      cmd_str = "cd " + pwd
      os.system(cmd_str)

  else:
    # or using the input scripts in the working directory
    input_list=['in.lj', 'in.rhodo', 'in.eam']
    iterate(input_list, config)

    

