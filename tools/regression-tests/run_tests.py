'''
UPDATE: Oct 25, 2023:
  Launching the LAMMPS binary under testing using a configuration defined in a yaml file (e.g. config.yaml)
  this way we can launch LAMMPS with mpirun with more flexibility. Also it simplifies the build configuration.

  python3 run_tests.py /path/to/lmp_binary

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
import sys
import re, yaml
import subprocess

try:
    from yaml import CSafeLoader as Loader
except ImportError:
    from yaml import SafeLoader as Loader

'''
  inputFileName:  input file with comments #REG:ADD and #REG:SUB as markers
  outputFileName: modified input file ready for testing
'''
def processing_markers(inputFileName, outputFileName):
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
  return the list of installed packages
'''
def get_installed_packages(lmp_binary):
  cmd_str = lmp_binary + " -h"
  p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
  output = p.stdout.split('\n')
  reading = False
  row = 0
  for l in output:
    if l != "":
      if l == "Installed packages:":
        reading = True
        n = row
      if reading == True and row > n:
        packages = l.strip()
        break
    row += 1
  return packages.split(" ")

'''
  launch LAMMPS using the configuration defined in the dictionary config with an input file
  TODO:
    - generate new reference values if needed
    - wrap subprocess with try/catch to handle exceptions
'''
def execute(lmp_binary, config, input_file_name, generate_ref_yaml=False):
  cmd_str = config['mpiexec'] + " " + config['mpiexec_numproc_flag'] + " " + config['nprocs'] + " "
  cmd_str += lmp_binary + " -in " + input_file_name + " " + config['args']
  print(f"Execute: {cmd_str}")
  p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
  output = p.stdout.split('\n')

'''
  Iterate over a list of input files
'''
def iterate(input_list):
  num_tests = len(input_list)
  test_id = 0
  # iterative over the input scripts
  for input in input_list:
    input_test=input + '.test'
    processing_markers(input, input_test)

    str_t = "\nRunning " + input_test + f" ({test_id+1}/{num_tests})"
    print(str_t)
    print(f"-"*len(str_t))

    # using the LAMMPS python module (for single-proc runs)
    #  lmp = lammps()
    #  lmp.file(input_test)
    
    # or more customizable with config.yaml
    execute(lmp_binary, config, input_test)

    # process thermo output
    thermo = extract_thermo("log.lammps")

    num_runs = len(thermo)
    if num_runs == 0:
       print(f"Failed with {input_test}\n")
       continue 

    print(f"Number of runs: {num_runs}")

    # read in the thermo yaml output
    thermo_ref_file = 'thermo.' + input + '.yaml'
    file_exist = os.path.isfile(thermo_ref_file)
    if file_exist == True:
        thermo_ref = extract_thermo(thermo_ref_file)
        # comparing output vs reference values
        width = 20
        print("Quantities".ljust(width) + "Output".center(width) + "Reference".center(width) + "Absolute Diff.".center(width))
        irun = 0
        num_fields = len(thermo[irun]['keywords'])
        for i in range(num_fields):
          val = thermo[irun]['data'][2][i]
          ref = thermo_ref[0]['data'][2][i]
          diff = abs(float(val) - float(ref))
          print(f"{thermo[0]['keywords'][i].ljust(width)} {str(val).rjust(20)} {str(ref).rjust(20)} {str(diff).rjust(20)}")
        print("-"*(4*width+3))
    else:
      print(f"{thermo_ref_file} does not exist")

    test_id = test_id + 1

'''
  TODO:
    - automate tagging the example input scripts of the installed packages
'''
if __name__ == "__main__":

  lmp_binary = ""
  # if lmp binary is specified in the command line
  if len(sys.argv) >= 2:
    lmp_binary = sys.argv[1]

  # read in the configuration of the tests
  with open("config.yaml", 'r') as f:
    config = yaml.load(f, Loader=Loader)
    print(f"Using configuration (config.yaml): {config}")
 
  # check if lmp_binary is specified in the config yaml
  if lmp_binary == "":
    if config['lmp_binary'] == "":
       print("Needs a valid LAMMPS binary")
       quit()
    else:
       lmp_binary = config['lmp_binary']

  packages = get_installed_packages(lmp_binary)
  print(f"List of installed packages: {packages}")

  # list of input scripts with markers #REG:SUB and #REG:ADD
  input_list=['in.lj', 'in.rhodo', 'in.eam']
  iterate(input_list)


  


