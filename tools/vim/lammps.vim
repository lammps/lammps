" Vim syntax file
" Language:	         Lammps Simulation Script File
" Maintainer:        Gerolf Ziegenhain <gerolf@ziegenhain.com>
" Updates:           Axel Kohlmeyer <akohlmey@gmail.com>
" Latest Revision:   2010-06-06

syn clear

syn keyword lammpsOutput      log write_restart dump undump thermo thermo_modify thermo_style print 
syn keyword lammpsRead        include read read_restart read_data
syn keyword lammpsLattice     boundary units atom_style lattice region create_box create_atoms dielectric
syn keyword lammpsLattice     delete_atoms change_box dimension newton replicate
syn keyword lammpsParticle    pair_coeff pair_style pair_modify mass angle_coeff angle_style atom_modify
syn keyword lammpsParticle    atom_style bond_coeff bond_style delete_bonds kspace_style kspace_modify
syn keyword lammpsParticle    dihedral_style dihedral_coeff improper_style improper_coeff
syn keyword lammpsSetup       min_style fix_modify run_style timestep neighbor fix unfix communicate
syn keyword lammpsSetup       neigh_modify reset_timestep velocity
syn keyword lammpsRun         minimize run  
syn keyword lammpsDefine      variable group

syn keyword lammpsRepeat      jump next loop

syn keyword lammpsOperator    equal add sub mult div 

syn keyword lammpsConditional if then else

syn region lammpsString			start=+'+ end=+'+	oneline
syn region lammpsString			start=+"+ end=+"+	oneline

syn match lammpsNumber		"\<[0-9]\+[ij]\=\>"
syn match lammpsFloat		"\<[0-9]\+\.[0-9]*\([edED][-+]\=[0-9]\+\)\=[ij]\=\>"
syn match lammpsFloat		"\.[0-9]\+\([edED][-+]\=[0-9]\+\)\=[ij]\=\>"
syn match lammpsFloat		"\<[0-9]\+[edED][-+]\=[0-9]\+[ij]\=\>"

syn match lammpsComment         "#.*$"

syn match lammpsVariable   "\$\({[a-zA-Z0-9]\+}\)"
syn match lammpsVariable   "\$[A-Za-z]"

if !exists("did_lammps_syntax_inits")
  let did_lammps_syntax_inits = 1
  hi link lammpsOutput        Function
  hi link lammpsRepeat        Repeat
  hi link lammpsRead          Include
  hi link lammpsLattice   		Typedef
  hi link lammpsParticle      Typedef
  hi link lammpsSetup         Typedef
  hi link lammpsDefine        Define
  hi link lammpsRun           Statement
  hi link lammpsNumber			Number
  hi link lammpsFloat			Float
  hi link lammpsString			String
  hi link lammpsComment			Comment
  hi link lammpsLoop          Repeat
  hi link lammpsVariable      Identifier
  hi link lammpsConditional   Conditional
  hi link lammpsOperator      Operator
endif

let b:current_syntax = "lammps"
