from pygments.lexer import RegexLexer, words, include, default
from pygments.token import *

LAMMPS_COMMANDS = ("angle_coeff", "angle_style", "atom_modify",
                   "atom_style", "angle_write", "balance", "bond_coeff",
                   "bond_style", "bond_write", "boundary", "clear",
                   "comm_modify", "comm_style", "compute_modify",
                   "create_atoms", "create_bonds", "create_box",
                   "delete_atoms", "delete_bonds", "dielectric",
                   "dihedral_coeff", "dihedral_style", "dihedral_write",
                   "dimension", "displace_atoms", "dump_modify",
                   "dynamical_matrix", "echo", "fitpod",
                   "fix_modify", "group2ndx", "hyper",
                   "improper_coeff", "improper_style", "include",
                   "info", "jump", "kim", "kspace_modify",
                   "kspace_style", "label", "labelmap", "lattice",
                   "log", "mass", "mdi", "message", "minimize",
                   "min_modify", "min_style", "molecule", "ndx2group",
                   "neb", "neb/spin", "neighbor", "neigh_modify",
                   "newton", "next", "package", "pair_coeff",
                   "pair_modify", "pair_style", "pair_write",
                   "partition", "plugin", "prd", "print", "processors",
                   "python", "quit", "read_data", "read_dump",
                   "read_restart", "replicate", "rerun", "reset_atoms",
                   "reset_timestep", "restart", "run", "run_style",
                   "server", "set", "shell", "special_bonds", "suffix",
                   "tad", "temper", "temper/grem", "temper/npt",
                   "thermo", "thermo_modify", "thermo_style",
                   "third_order", "timer", "timestep", "units",
                   "velocity", "write_coeff", "write_data",
                   "write_restart")

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
            (r'fix_modify\s+', Keyword, 'modify_cmd'),
            (r'compute\s+', Keyword, 'compute'),
            (r'compute_modify\s+', Keyword, 'modify_cmd'),
            (r'dump\s+', Keyword, 'dump'),
            (r'dump_modify\s+', Keyword, 'modify_cmd'),
            (r'region\s+', Keyword, 'region'),
            (r'variable\s+', Keyword, 'variable_cmd'),
            (r'group\s+', Keyword, 'group'),
            (r'change_box\s+', Keyword, 'change_box'),
            (r'create_box\s+', Keyword, 'create_box'),
            (r'delete_bonds\s+', Keyword, 'id_cmd'),
            (r'displace_atoms\s+', Keyword, 'id_cmd'),
            (r'dynamical_matrix\s+', Keyword, 'id_cmd'),
            (r'group2ndx\s+', Keyword, 'ndx_cmd'),
            (r'ndx2group\s+', Keyword, 'ndx_cmd'),
            (r'jump\s+', Keyword, 'jump_cmd'),
            (r'label\s+', Keyword, 'jump_cmd'),
            (r'next\s+', Keyword, 'id_cmd'),
            (r'kim\s+', Keyword, 'kim_cmd'),
            (r'uncompute\s+', Keyword, 'id_cmd'),
            (r'unfix\s+', Keyword, 'id_cmd'),
            (r'undump\s+', Keyword, 'id_cmd'),
            (r'velocity\s+', Keyword, 'id_cmd'),
            (r'write_coeff\s+', Keyword, 'ndx_cmd'),
            (r'write_data\s+', Keyword, 'ndx_cmd'),
            (r'write_dump\s+', Keyword, 'write_dump'),
            (r'write_restart\s+', Keyword, 'ndx_cmd'),
            include('conditionals'),
            include('keywords'),
            (r'#.*?\n', Comment),
            (r'"', String, 'string'),
            (r'\'', String, 'single_quote_string'),
            (r'[0-9]+:[0-9]+(:[0-9]+)?', Number),
            (r'[0-9]+(\.[0-9]+)?([eE]\-?[0-9]+)?', Number),
            (r'\$?\(', Name.Variable, 'expression'),
            (r'\$\{', Name.Variable, 'variable'),
            (r'[\w_\.\[\]]+', Name),
            (r'\$[\w_]+', Name.Variable),
            (r'\s+', Whitespace),
            (r'[\+\-\*\^\|\/\!%&=<>]', Operator),
            (r'[\~\.\w_:,@\-\/\\0-9]+', Text),
        ],
        'conditionals' : [
            (words(('if','else','elif','then'), suffix=r'\b', prefix=r'\b'), Keyword)
        ]
        ,
        'keywords' : [
            (words(LAMMPS_COMMANDS, suffix=r'\b', prefix=r'^\s*'), Keyword)
        ]
        ,
        'variable' : [
            (r'[^\}]+', Name.Variable),
            (r'\}', Name.Variable, '#pop'),
        ],
        'string' : [
            (r'[^"]+', String),
            (r'"', String, '#pop'),
        ],
        'single_quote_string' : [
            (r'[^\']+', String),
            (r'\'', String, '#pop'),
        ],
        'expression' : [
            (r'[^\(\)]+', Name.Variable),
            (r'\(', Name.Variable, 'expression'),
            (r'\)', Name.Variable, '#pop'),
        ],
        'modify_cmd' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
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
        'variable_cmd' : [
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
        'create_box' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            (r'\s+', Whitespace, 'group_id'),
            default('#pop')
        ],
        'id_cmd' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'ndx_cmd' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            default('#pop')
        ],
        'jump_cmd' : [
            (r'[\w_\-\.\[\]]+', Literal.String.Char),
            default('#pop')
        ],
        'kim_cmd' : [
            (r'[\w_\-\.\[\]]+', Literal.String.Single),
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
