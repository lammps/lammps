## reorder_remd_traj

LAMMPS Replica Exchange Molecular Dynamics (REMD) trajectories (implemented using the temper command) are arranged by replica, i.e., each trajectory is a continuous replica that records all the ups and downs in temperature. However, often the requirement is  that trajectories be continuous in temperature. This requires the LAMMPS REMD trajectories to be re-ordered, which LAMMPS does not do automatically. (see the discussion [here](https://lammps.sandia.gov/threads/msg60440.html)). The reorderLAMMPSREMD tool does exactly this in parallel (using MPI)

(Protein folding trajectories in [Sanyal, Mittal and Shell, JPC, 2019, 151(4), 044111](https://aip.scitation.org/doi/abs/10.1063/1.5108761) were ordered in temperature space using this tool)

#### Author

Tanmoy Sanyal, Shell lab, UC Santa Barbara

(currently at UC San Francisco)

email: tanmoy dot 7989 at gmail.com

#### Features

- reorder LAMMPS REMD trajectories by temperature keeping only desired frames.
  Note: this only handles LAMMPS format trajectories (i.e., lammpstrj format)
  Trajectories can be gzipped or bz2-compressed. The trajectories are assumed to
  be named as \<prefix>\.%d.lammpstrj[.gz or .bz2]

- (optionally) calculate configurational weights for each frame at each
  temperature if potential energies are supplied (only implemented for the canonical (NVT) ensemble)

#### Dependencies

[`mpi4py`](https://mpi4py.readthedocs.io/en/stable/)  
[`pymbar`](https://pymbar.readthedocs.io/en/master/) (for getting configurational weights)  
[`tqdm`](https://github.com/tqdm/tqdm) (for printing pretty progress bars)  
[`StringIO`](https://docs.python.org/2/library/stringio.html) (or [`io`](https://docs.python.org/3/library/io.html) if in Python 3.x)

#### Example

###### REMD Simulation specs 
Suppose you ran a REMD simulation for the peptide example using the CHARMM forcefield (see lammps/examples/peptide) in Lammps with the following settings:

- number of replicas = 16
- temperatures used (in K): 200 209 219 230 241 252 264 276 289 303 317 332 348 365 382 400 (i.e., exponentially distributed in the range 270-400 K)
- timestep = 2 fs
- total number of timesteps simulated using temper = 2000 (i.e. 4 ps)
- swap frequency = temperatures swapped after every this many steps = `ns` = 10 (i.e. 20 fs)
- write frequency = trajectory frame written to disk after this many steps (using the dump command) = `nw` = 20 (i.e. 40 fs)

###### LAMMPS output
So, when the dust settles,

- You'll have 16 replica trajectories. For this tool to work, each replica traj must be named: `<prefix>.<n>.lammpstrj[.gz or .bz2]`, where,
  - `prefix` = some common prefix for all your trajectories and (say it is called "peptide")` 
  - `n` = replica number (0-15 in this case). Note: trajectories **must be in default LAMMPS format **(so stuff like dcd won't work)

- You will also have a master LAMMPS log file (`logfn`) that contains the swap history of all the replicas
  (for more details see [here](https://lammps.sandia.gov/doc/temper.html). Assume that this is called `log.peptide`

- Further you must have a txt file that numpy can read which stores all the temperature values (say this is called `temps.txt`)

######  Your desired output
- The total number of timesteps you want consider as production (i.e. after equilbration)  = 1000 (i.e. last 2 ps)

- Reordered trajectories at temperatures 200 K, 276 K, 400 K

- Configurational log-weight calculation (using [`pymbar`](https://github.com/choderalab/pymbar)). Here, this is limited to the canonical (NVT) ensemble **and without biasing restraints** in your simulation. To do this you'd need to have a file (say called `ene.dat`) that stores a 2D  (K X N) array of total potential energies, where,

  - K = total number of replicas = 16, and N = total number of frames in each replica trajectory (= 1000 / 20 = 50 in this case) 

  - `ene[k,n]` = energy from n-th frame of k-th replica.

###### Using the tool (description of the workflow)
Assume you have 16 processors at your disposal. When you run the following:

```bash
mpirun -np 16 python reorder_remd_traj.py peptide -logfn log.peptide -tfn temps.txt -ns 10 -nw 20 -np 1000 -ot 200 276 400 -logw -e ene.peptide -od ./output
```

1. First the temperature swap history file (`log.peptide` in this case) is read. This is done on one processor since it is usually fast.
2. Then the (compressed or otherwise) LAMMPS replica trajectories are read in parallel. So if you have less processors than replicas at this stage, it'll be slower.
3. Then using the frame ordering generated in (1), trajectory frames read in (2) are re-ordered and written to disk in parallel. Each processor writes one trajectory. So, If you request reordered trajectories for less temperatures (3 in this case) than the total number of temperatures (16), then 16-3 = 13 processors will be retired.
4. If you have further requested configurational log-weight calculation, then they will be done on a single processor since pymbar is pretty fast.
5. Finally you will have 3 LAMMPS trajectories of the form ``peptide.<temp>.lammpstrj.gz`` each with 1000 / 20 = 50 frames,  where `<temp>` = 200, 276, 400. If you request reordering at a temperature like say 280 K which is not present in the supplied temp schedule (as written in `temps.txt`), the closest temperature (276 K) will be chosen.

For more details, use the help menu generated by the tool by using:
python reorder_remd_traj.py -h

###### Caveats
- This tool crawls through the replica trajectories and creates index files that are hidden. These are called .byteind_`<replicanum>'.gz files. You may delete these if you want, but subsequent replica reads will be slow in that case.

- When writing trajectories to disk, the trajectories are first written to a buffer in memory, and then finally dumped all-at-once to the disk. While this makes the tool very fast, it can cause out-of-memory errors for very large trajectories. A useful feature might be to write to the buffer in batches and emptying to disk when some (predefined) max-buffer-size is exceeded.
