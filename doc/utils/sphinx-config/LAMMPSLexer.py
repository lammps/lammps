from pygments.lexer import RegexLexer, words, include, default
from pygments.token import *

LAMMPS_COMMANDS = ("angle_coeff", "angle_style", "atom_modify", "atom_style",
"balance", "bond_coeff", "bond_style", "bond_write", "boundary", "box",
"clear", "comm_modify", "comm_style",
"compute_modify", "create_atoms", "create_bonds", "create_box", "delete_atoms",
"delete_bonds", "dielectric", "dihedral_coeff", "dihedral_style", "dimension",
"displace_atoms", "dump_modify", "dynamical_matrix", "echo",
"fix_modify", "group2ndx", "hyper", "if", "improper_coeff",
"improper_style", "include", "info", "jump", "kim_init", "kim_interactions",
"kim_param", "kim_query", "kspace_modify", "kspace_style", "label", "lattice",
"log", "mass", "message", "minimize", "min_modify", "min_style", "molecule",
"ndx2group", "neb", "neb/spin", "neighbor", "neigh_modify", "newton", "next",
"package", "pair_coeff", "pair_modify", "pair_style", "pair_write",
"partition", "prd", "print", "processors", "python", "quit", "read_data",
"read_dump", "read_restart", "replicate", "rerun", "reset_ids",
"reset_timestep", "restart", "run", "run_style", "server", "set", "shell",
"special_bonds", "suffix", "tad", "temper", "temper/grem", "temper/npt",
"thermo", "thermo_modify", "thermo_style", "then", "third_order", "timer", "timestep",
"units", "velocity", "write_coeff",
"write_data", "write_restart")

#fix ID group-ID style args
#compute ID group-ID style args
#dump ID group-ID style N file args
#region ID style args keyword arg ...
#variable name style args ...
#group ID style args
#uncompute compute-ID
#undump dump-ID
#unfix fix-ID
#write_dump group-ID style file dump-args modify dump_modify-args

class LAMMPSLexer(RegexLexer):
    name = 'LAMMPS'
    tokens = {
        'root': [
            (r'fix\s+', Keyword, 'fix'),
            (r'compute\s+', Keyword, 'compute'),
            (r'dump\s+', Keyword, 'dump'),
            (r'region\s+', Keyword, 'region'),
            (r'variable\s+', Keyword, 'variable'),
            (r'group\s+', Keyword, 'group'),
            (r'change_box\s+', Keyword, 'change_box'),
            (r'uncompute\s+', Keyword, 'uncompute'),
            (r'unfix\s+', Keyword, 'unfix'),
            (r'undump\s+', Keyword, 'undump'),
            (r'write_dump\s+', Keyword, 'write_dump'),
            include('keywords'),
            (r'#.*?\n', Comment),
            ('"', String, 'string'),
            ('\'', String, 'single_quote_string'),
            (r'[0-9]+(\.[0-9]+)?([eE]\-?[0-9]+)?', Number),
            ('\$?\(', Name.Variable, 'expression'),
            ('\$\{', Name.Variable, 'variable'),
            (r'[\w_\.\[\]]+', Name),
            (r'\$[\w_]+', Name.Variable),
            (r'\s+', Whitespace),
            (r'[\+\-\*\^\|\/\!%&=<>]', Operator),
        ],
        'keywords' : [
            (words(LAMMPS_COMMANDS, suffix=r'\b', prefix=r'^'), Keyword)
        ]
        ,
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
        ],
        'fix' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            (r'\s+', Whitespace, 'group_id'),
            default('#pop')
        ],
        'compute' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            (r'\s+', Whitespace, 'group_id'),
            default('#pop')
        ],
        'dump' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            (r'\s+', Whitespace, 'group_id'),
            default('#pop')
        ],
        'region' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'variable' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'group' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'change_box' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'unfix' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'undump' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'uncompute' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'write_dump' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'group_id' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop:2')
        ]
    }
