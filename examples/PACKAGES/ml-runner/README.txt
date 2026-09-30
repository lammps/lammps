This directory contains an example simulation for the ML-RUNNER package. 
It runs a short NVT molecular dynamics simulation of 192 bulk water molecules.

Files included:
in.ml-runner.H2O             # LAMMPS input script for the NVT simulation
H2O.data                     # LAMMPS data file containing 64 bulk water molecules

2G-H2O-HDNNP/                # Directory containing the 2G-HDNNP trained on bulk water reference data
  ├── input.nn               # RuNNer settings, network architecture and feature map definition
  ├── scaling.data           # RuNNer feature map scaling information
  ├── weights_short.001.data # RuNNer weights for the Hydrogen atomic neural network
  └── weights_short.008.data # RuNNer weights for the Oxygen atomic neural network

To run this example, you must compile LAMMPS with the ML-RUNNER package enabled
and the RuNNer library linked. Then, execute the following command:

lmp -in in.ml-runner.H2O
