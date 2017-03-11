#!/usr/bin/env python
"""
   Unfortunately, the lines in a PDB files are not always listed in the 
   correct order.  Software that reads PDB files is expected to re-sort this
   data before interpreting it.  (One reason I don't like them.)
   This script reads a PDB file from the standard input, sorts the lines
   according to ChainID, SeqNum (residue-number), Icode (insert-code), 
   and AtomID (in that order) and prints the result to the standard-out.
   Only the ATOM and HETATM records are effected.
   All other lines in the PDB file are printed back to the user verbatim.
   Note: the "altLoc" column (character #17 on each line) is ignored
         and is not used for sorting purposes.

"""

from collections import defaultdict


# In order to specify any amino acid in a PDB file, you must provide 3
# identifiers:
#   the ChainID a single letter specifying
#   the SeqNum an integer indicating the location within that chain
#   the ICode ("insert code" usually 0.  I don't know why this number is 
#              necessary. ..For the record, I never loved the PDB file format.)


class AtomDescr:
    def __init__(self, setChainID, setSeqNum, setICode, setAtomID):
        self.chainID = setChainID
        self.seqNum = setSeqNum
        self.iCode = setICode
        self.atomID = setAtomID

    #I plan to store this information in a python dictionary.
    #Unfortunately, in order to use such classes as keys in python dictionaries
    #I must define comparison operators, and a hash function.
    #In retrospect I figured out it would have been easier just to use tuples
    #as dictionary keys.  I suppose it was a good excercise to try to do it
    #using python classes instead.  I'm sorry it made the code so long.

    def __le__(self, other):
        #return ((self.chainID < other.chainID) or ((self.chainID == other.chainID) and ((self.seqNum < other.seqNum) or ((self.seqNum == other.seqNum) and (self.iCode <= other.iCode)))))))
        # instead I'll exploit python's ability to compare tuples
        return (self.chainID, self.seqNum, self.iCode, self.atomID) <= (other.chainID, other.seqNum, other.iCode, other.atomID)

    def __lt__(self, other):
        return (self.chainID, self.seqNum, self.iCode, self.atomID) < (other.chainID, other.seqNum, other.iCode, other.atomID)

    def __eq__(self, other):
        #return ((self.chainID == x.chainID) and (self.seqNum == x.seqNum) and (self.iCode == x.iCode))
        return (self.chainID, self.seqNum, self.iCode, self.atomID) == (other.chainID, other.seqNum, other.iCode, other.atomID)

    def __ne__(self, other):
        return not __eq__(self, other)

    def __gt__(self, other):
        return not __le__(self, other)

    def __ge__(self, other):
        return not __lt__(self, other)

    def __cmp__(self, other):
        if __lt__(self, other):
            return -1
        elif __gt__(self, other):
            return 1
        else:
            return 0

    def __hash__(self):
        numChainIDs = 128
        numICodes = 128
        i = self.seqNum
        i *= numChainIDs
        i += ord(self.chainID)
        i *= numICodes
        i += ord(self.iCode)
        i *=10
        i += self.atomID
        return i




import sys
from operator import attrgetter
g_program_name = __file__.split('/')[-1]
g_version_str = 0.11
g_date_str = 2013-9-18

if len(sys.argv) == 1:
    use_all_residues = True
elif len(sys.argv) == 7:
    use_all_residues = False
    first = AtomDescr(sys.argv[1], int(sys.argv[2]), sys.argv[3], 0)
    last  = AtomDescr(sys.argv[4], int(sys.argv[5]), sys.argv[6], 2147483647)
else:
    sys.stderr.write("Error("+g_program_name+"): This program requires either 0 or 6 arguments.\n"
                     "       By default, the the sequence is extracted from the entire PDB file.\n"
                     "       In that case, no arguments are required.\n"
                     "       Alternately, you can limit the selection to a single interval of\n"
                     "       residues from one of the chains in the PDB file.\n"
                     "       To specify an interval, you must passing 6 arguments to this program.\n"
                     "       This program requires a pair of residues to designate the first and\n"
                     "       last members of the interval.  Each residue requires 3 identifiers.\n"
                     "       Consequently the six arguments needed are:\n"
                     "ChainID_first SeqNum_first ICode_first ChainID_last SeqNum_last ICode_last\n")
    exit(-1)


atoms2lines = defaultdict(list)

for line in sys.stdin:
    if (line[0:6] == "ATOM  ") or (line[0:6] == "HETATM"):
        atomID    = int(line[6:11])
        #atomType  = line[12:16]
        #altLoc    = line[16:17]
        iCode     = line[26:27]
        #resType   = line[17:20]
        chainID   = line[21:22]
        seqNumStr = line[22:26]
        seqNum    = int(seqNumStr)
        atomdescr = AtomDescr(chainID, int(seqNumStr), iCode, int(atomID))
        atoms2lines[atomdescr].append(line.rstrip('\n'))
    else:
        sys.stdout.write(line)

# Extract an (unordered) list of the atomdescrs of the atoms in the sequence
atomdescrs = [atomdescr for atomdescr in atoms2lines]

# Residues in PDB files are often not listed in order.
# Consequently, we must sort the list by chainID, seqNum, and finnaly iCode:
sequence_of_atomdescrs = sorted(atomdescrs, key=attrgetter('chainID','seqNum','iCode','atomID'))

for atomdescr in sequence_of_atomdescrs:
    for line in atoms2lines[atomdescr]:
        sys.stdout.write(line+'\n')

