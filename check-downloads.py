#!/usr/bin/env python3

import subprocess
import json

try:
    output = subprocess.run('gh api https://api.github.com/repos/lammps/lammps/releases',
                   shell=True, capture_output=True)
except subprocess.CalledProcessError as e:
    print("API call failed with:", e.output)

releases = json.loads(output.stdout)
print("Recent releases: ", len(releases))
for rel in releases:
    if len(rel['assets']) > 0:
        print("Release: ", rel['name'])
        for asset in rel['assets']:
            print("Asset: %s  Downloads: %d" % (asset['name'], asset['download_count']))
