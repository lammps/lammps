from pygments.lexer import RegexLexer, words
from pygments.token import *

LAMMPS_COMMANDS = ("angle_coeff", "angle_style", "atom_modify", "atom_style",
"balance", "bond_coeff", "bond_style", "bond_write", "boundary", "box",
"change_box", "clear", "comm_modify", "comm_style", "compute",
"compute_modify", "create_atoms", "create_bonds", "create_box", "delete_atoms",
"delete_bonds", "dielectric", "dihedral_coeff", "dihedral_style", "dimension",
"displace_atoms", "dump", "dump_modify", "dynamical_matrix", "echo", "fix",
"fix_modify", "group", "group2ndx", "hyper", "if", "improper_coeff",
"improper_style", "include", "info", "jump", "kim_init", "kim_interactions",
"kim_param", "kim_query", "kspace_modify", "kspace_style", "label", "lattice",
"log", "mass", "message", "minimize", "min_modify", "min_style", "molecule",
"ndx2group", "neb", "neb/spin", "neighbor", "neigh_modify", "newton", "next",
"package", "pair_coeff", "pair_modify", "pair_style", "pair_write",
"partition", "prd", "print", "processors", "python", "quit", "read_data",
"read_dump", "read_restart", "region", "replicate", "rerun", "reset_ids",
"reset_timestep", "restart", "run", "run_style", "server", "set", "shell",
"special_bonds", "suffix", "tad", "temper", "temper/grem", "temper/npt",
"thermo", "thermo_modify", "thermo_style", "then", "third_order", "timer", "timestep",
"uncompute", "undump", "unfix", "units", "variable", "velocity", "write_coeff",
"write_data", "write_dump", "write_restart")

class LAMMPSLexer(RegexLexer):
    name = 'LAMMPS'
    tokens = {
        'root': [
            (words(LAMMPS_COMMANDS, suffix=r'\b', prefix=r'^'), Keyword),
            (r'#.*?\n', Comment),
            ('"', String, 'string'),
            ('\'', String, 'single_quote_string'),
            (r'[0-9]+(\.[0-9]+)?([eE]\-?[0-9]+)?', Number),
            ('\$?\(', Name.Variable, 'expression'),
            ('\$\{', Name.Variable, 'variable'),
            (r'[\w_\.\[\]]+', Name),
            (r'\$[\w_]+', Name.Variable),
            (r'\s+', Whitespace),
            (r'[\+\-\*\/&=<>]', Operator),
        ],
        'variable' : [
            ('[^\}]+', Name.Variable),
            ('\}', Name.Variable, '#pop'),
        ],
        'string' : [
            ('[^"]+', String),
            ('"', String, '#pop'),
        ],
        'single_quote_string' : [
            ('[^\']+', String),
            ('\'', String, '#pop'),
        ],
        'expression' : [
            ('[^\(\)]+', Name.Variable),
            ('\(', Name.Variable, 'expression'),
            ('\)', Name.Variable, '#pop'),
        ]
    }
