" Vim syntax file
" Language:          Lammps Simulation Script File
" Maintainer:        Gerolf Ziegenhain <gerolf@ziegenhain.com>
" Updates:           Axel Kohlmeyer <akohlmey@gmail.com>, Sam Bateman <sam.bateman@nrlssc.navy.mil>, Daniel Möller Montull <d.moller.m@gmail.com>, Eryk Skalinski <eskalinski@protonmail.com>
" Latest Revision:   2023-07-15

syn clear

" Add '/' to list of valid keyword characters
set iskeyword+=/

syn keyword     lammpsOutput    log write_data write_dump write_coeff info shell write_restart restart dump undump thermo thermo_modify
syn keyword     lammpsOutput    thermo_style print timer
syn keyword     lammpsRead      include read_restart read_data read_dump molecule
syn keyword     lammpsLattice   boundary units atom_style lattice region create_box create_atoms dielectric
syn keyword     lammpsLattice   delete_atoms displace_atoms change_box dimension replicate
syn keyword     lammpsParticle  pair_coeff pair_style pair_modify pair_write mass velocity angle_coeff angle_style angle_write
syn keyword     lammpsParticle  atom_modify atom_style bond_coeff bond_style bond_write create_bonds delete_bonds kspace_style
syn keyword     lammpsParticle  kspace_modify dihedral_style dihedral_coeff dihedral_write improper_style improper_coeff labelmap
syn keyword     lammpsSetup     min_style min_modify fix_modify run_style timestep neighbor neigh_modify fix unfix suffix special_bonds dump_modify
syn keyword     lammpsSetup     balance box clear comm_modify comm_style newton package processors reset_atoms reset_ids reset_timestep
syn keyword     lammpsRun       minimize minimize/kk run rerun tad neb neb/spin prd quit server temper/npt temper/grem temper
syn keyword     lammpsRun       message hyper dynamical_matrix dynamical_matrix/kk third_order third_order/kk fitpod
syn keyword     lammpsDefine    variable group compute python set uncompute kim_query kim group2ndx ndx2group mdi

syn keyword     lammpsRepeat    jump next loop label

syn keyword     lammpsOperator  equal add sub mult div

syn keyword     lammpsConditional if then elif else

syn keyword     lammpsSpecial   EDGE NULL INF &

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
