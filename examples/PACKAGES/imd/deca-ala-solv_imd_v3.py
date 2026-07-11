"""
For use with 'in.deca-ala-solv_imd_v3'.

Tested with imdclient v0.1.4 and MDAnalysis v2.8.0
"""
from imdclient.IMD import IMDReader
import MDAnalysis as mda

u = mda.Universe('data.deca-ala-solv', "imd://localhost:5678", topology_format='DATA')

for ts in u.trajectory:
    print(ts.time)
    print(ts.velocities)