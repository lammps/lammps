" Vim syntax file
" Language:              Lammps Simulation Script File
" Maintainer:        Gerolf Ziegenhain <gerolf@ziegenhain.com>
" Updates:           Axel Kohlmeyer <akohlmey@gmail.com>, Sam Bateman <sam.bateman@nrlssc.navy.mil>, Daniel Möller Montull <d.moller.m@gmail.com>, Eryk Skalinski <eskalinski@protonmail.com>
" Latest Revision:   2022-08-17

syn clear

syn keyword     lammpsOutput    log write_data write_dump info shell write_restart restart dump undump thermo thermo_modify thermo_style print timer
syn keyword     lammpsRead      include read_restart read_data read_dump molecule
syn keyword     lammpsLattice   boundary units atom_style lattice region create_box create_atoms dielectric
syn keyword     lammpsLattice   delete_atoms displace_atoms change_box dimension replicate
syn keyword     lammpsParticle  pair_coeff pair_style pair_modify pair_write mass velocity angle_coeff angle_style
syn keyword     lammpsParticle  atom_modify atom_style bond_coeff bond_style bond_write create_bonds delete_bonds kspace_style
syn keyword     lammpsParticle  kspace_modify dihedral_style dihedral_coeff improper_style improper_coeff
syn keyword     lammpsSetup     min_style fix_modify run_style timestep neighbor neigh_modify fix unfix suffix special_bonds
syn keyword     lammpsSetup     balance box clear comm_modify comm_style newton package processors reset_ids reset_timestep
syn keyword     lammpsRun       minimize run rerun tad neb prd quit server temper temper/grem temper/npt
syn keyword     lammpsRun       min/spin message hyper dynamical_matrix
syn keyword     lammpsDefine    variable group compute python set uncompute kim_query

syn keyword     lammpsRepeat    jump next loop

syn keyword     lammpsOperator  equal add sub mult div

syn keyword     lammpsConditional if then elif else

syn keyword     lammpsSpecial   EDGE NULL &

syn region      lammpsString    start=+'+ end=+'+       oneline
syn region      lammpsString    start=+"+ end=+"+       oneline

syn match       lammpsNumber    "\<[0-9]\+[ij]\=\>"
syn match       lammpsFloat     "\<[0-9]\+\.[0-9]*\([edED][-+]\=[0-9]\+\)\=[ij]\=\>"
syn match       lammpsFloat     "\.[0-9]\+\([edED][-+]\=[0-9]\+\)\=[ij]\=\>"
syn match       lammpsFloat     "\<[0-9]\+[edED][-+]\=[0-9]\+[ij]\=\>"

syn match       lammpsComment   "#\(.*&\s*\n\)*.*$"

syn match       lammpsVariable  "\$\({[a-zA-Z0-9_]\+}\)"
syn match       lammpsVariable  "\$[A-Za-z]"

if !exists("did_lammps_syntax_inits")
  let did_lammps_syntax_inits = 1
  hi def link lammpsOutput              Function
  hi def link lammpsRepeat              Repeat
  hi def link lammpsRead                Include
  hi def link lammpsLattice             Typedef
  hi def link lammpsParticle            Typedef
  hi def link lammpsSetup               Typedef
  hi def link lammpsDefine              Define
  hi def link lammpsRun                 Statement
  hi def link lammpsNumber              Number
  hi def link lammpsFloat               Float
  hi def link lammpsString              String
  hi def link lammpsComment             Comment
  hi def link lammpsLoop                Repeat
  hi def link lammpsVariable            Identifier
  hi def link lammpsConditional         Conditional
  hi def link lammpsOperator            Operator
  hi def link lammpsSpecial             Number
endif

let b:current_syntax = "lammps"
" vim: set expandtab
