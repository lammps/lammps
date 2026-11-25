from pygments.lexer import RegexLexer, words, include, default, bygroups
from pygments.token import *

LAMMPS_COMMANDS = ("angle_coeff", "angle_style", "angle_write", "atom_modify",
                   "atom_style", "balance", "bond_coeff", "bond_style", "bond_write",
                   "boundary", "change_box", "clear", "comm_modify", "comm_style",
                   "compute_modify", "create_atoms", "create_bonds", "create_box",
                   "delete_atoms", "delete_bonds", "dielectric", "dihedral_coeff",
                   "dihedral_style", "dihedral_write", "dimension", "displace_atoms",
                   "dump_modify", "dynamical_matrix", "echo", "fenix", "fitpod",
                   "fix_modify", "geturl", "group2ndx", "hyper", "improper_coeff",
                   "improper_style", "include", "info", "jump", "kim", "kspace_modify",
                   "kspace_style", "label", "labelmap", "lattice", "log", "mass",
                   "mdi", "minimize", "min_modify", "min_style", "molecule",
                   "ndx2group", "neb", "neb/spin", "neighbor", "neigh_modify",
                   "newton", "next", "package", "pair_coeff", "pair_modify",
                   "pair_style", "pair_write", "partition", "plugin", "prd", "print",
                   "processors", "python", "quit", "read_data", "read_dump",
                   "read_restart", "region2vmd", "replicate", "rerun", "reset_atoms",
                   "reset_timestep", "restart", "run", "run_style", "set", "shell",
                   "special_bonds", "suffix", "tad", "temper", "temper/grem",
                   "temper/npt", "thermo", "thermo_modify", "thermo_style",
                   "third_order", "timer", "timestep", "units", "velocity",
                   "write_coeff", "write_data", "write_molecule", "write_restart")

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

# rule for a command at the beginning of a line whose following word(s)
# are handled by a dedicated state (IDs, group-IDs, labels, ...); anchoring
# at the line start avoids false matches on keyword arguments of the same
# name (e.g. "region" or "dump" as arguments of other commands)
def cmd_rule(command, state):
    return (r'(^[ \t]*)(' + command + r')(\s+)',
            bygroups(Whitespace, Keyword, Whitespace), state)

class LAMMPSLexer(RegexLexer):
    name = 'LAMMPS'
    tokens = {
        'root': [
            cmd_rule('fix', 'define_cmd'),
            cmd_rule('fix_modify', 'modify_cmd'),
            cmd_rule('compute', 'define_cmd'),
            cmd_rule('compute_modify', 'modify_cmd'),
            cmd_rule('dump', 'define_cmd'),
            cmd_rule('dump_modify', 'modify_cmd'),
            cmd_rule('region', 'region'),
            cmd_rule('variable', 'variable_cmd'),
            cmd_rule('group', 'group'),
            cmd_rule('change_box', 'change_box'),
            cmd_rule('create_box', 'create_box'),
            cmd_rule('delete_bonds', 'id_cmd'),
            cmd_rule('displace_atoms', 'id_cmd'),
            cmd_rule('dynamical_matrix', 'id_cmd'),
            cmd_rule('group2ndx', 'ndx_cmd'),
            cmd_rule('ndx2group', 'ndx_cmd'),
            cmd_rule('jump', 'jump_cmd'),
            cmd_rule('label', 'jump_cmd'),
            cmd_rule('next', 'id_cmd'),
            cmd_rule('kim', 'kim_cmd'),
            cmd_rule('uncompute', 'id_cmd'),
            cmd_rule('unfix', 'id_cmd'),
            cmd_rule('undump', 'id_cmd'),
            cmd_rule('velocity', 'id_cmd'),
            cmd_rule('write_coeff', 'ndx_cmd'),
            cmd_rule('write_data', 'ndx_cmd'),
            cmd_rule('write_dump', 'write_dump'),
            cmd_rule('write_restart', 'ndx_cmd'),
            cmd_rule('angle_style', 'style_name'),
            cmd_rule('atom_style', 'style_name'),
            cmd_rule('bond_style', 'style_name'),
            cmd_rule('dihedral_style', 'style_name'),
            cmd_rule('improper_style', 'style_name'),
            cmd_rule('kspace_style', 'style_name'),
            cmd_rule('min_style', 'style_name'),
            cmd_rule('pair_style', 'style_name'),
            cmd_rule('run_style', 'style_name'),
            include('conditionals'),
            include('keywords'),
            (r'#.*?\n', Comment),
            (r'&[ \t]*\n', Literal.String.Char),
            (r'"', String, 'string'),
            (r'\'', String, 'single_quote_string'),
            # special argument words accepted by several commands
            (words(('NULL', 'EDGE', 'INF', 'SELF'), prefix=r'\b', suffix=r'\b'),
             Keyword.Constant),
            # references to variables (v_), computes (c_), fixes (f_), and
            # custom per-atom properties (i_, d_, i2_, d2_)
            (r'\b(i2|d2|[cdfiv])_[\w\[\]]+', Name.Variable),
            (r'[0-9]+:[0-9]+(:[0-9]+)?', Number),
            (r'([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?', Number),
            (r'\$?\(', Name.Variable, 'expression'),
            (r'\$\{', Name.Variable, 'variable'),
            # words with forward slashes (style names, file paths) are one name
            (r'[A-Za-z][\w\.\[\]]*(/[\w\.\[\]]+)+', Name),
            (r'[\w_\.\[\]]+', Name),
            (r'\$[\w_]+', Name.Variable),
            (r'[^\S\n]+', Whitespace),
            (r'\n', Whitespace),
            (r'[\+\-\*\^\|\/\!%&=<>]', Operator),
            (r'[\~\.\w_:,@\-\/\\0-9]+', Text),
        ],
        'conditionals' : [
            (words(('if','else','elif','then'), suffix=r'\b', prefix=r'\b'), Keyword)
        ]
        ,
        'keywords' : [
            (words(LAMMPS_COMMANDS, suffix=r'\b', prefix=r'^[ \t]*'), Keyword)
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
        # fix/compute/dump: ID, then group-ID, then style name
        'define_cmd' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            (r'[^\S\n]+', Whitespace, ('#pop', 'style_name', 'group_id')),
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
            (r'[^\S\n]+', Whitespace, ('#pop', 'id_cmd')),
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
        # write_dump: group-ID, then style name
        'write_dump' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            (r'[^\S\n]+', Whitespace, ('#pop', 'style_name')),
            default('#pop')
        ],
        'group_id' : [
            (r'[\w_\-\.\[\]]+', Name.Variable.Identifier),
            (r'[^\S\n]+', Whitespace, '#pop'),
            default('#pop:2')
        ],
        # one style name; entered with the style word as the next token
        'style_name' : [
            (r'[A-Za-z][\w/]*', Name.Class, '#pop'),
            default('#pop')
        ]
    }
