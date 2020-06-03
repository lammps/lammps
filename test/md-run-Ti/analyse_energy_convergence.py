import math


def mean(lst):
    return sum(lst) / len(lst)


def std(lst):
    m = mean(lst)
    var = mean([(el - m) ** 2 for el in lst])
    return math.sqrt(var)


import sys

with open(sys.argv[1]) as f:
    lines = f.readlines()

found = False
for i, line in enumerate(lines):
    if 'Step TotEng PotEng KinEng Pxx Pyy Pzz Press Temp Volume' in line:
        found = True
        break

i += 1

if not found:
    raise ValueError("Begin of MD trajectory not found")

found = False
for j in range(i, len(lines)):
    line = lines[j]
    if 'Loop time of' in line:
        found = True
        break

if not found:
    raise ValueError("End of MD trajectory not found")

data_lines = lines[i:j]

data_lines = list(map(lambda s: s.strip().split(), data_lines))

data = [[float(e) for e in line] for line in data_lines]
try:
    n_atoms = int(lines[j].strip().split()[-2])
    print("Number of atoms: ", n_atoms)
except Exception as e:
    print("Error: ", e)
en_data = [l[1] for l in data]
std_enrgy_deviation = std(en_data) / abs(mean(en_data))

print("Std. of NVE energy deviation: ", std_enrgy_deviation)
print("Std. of NVE energy deviation per atom: ", std_enrgy_deviation / n_atoms)
if std_enrgy_deviation > 4e-07:
    # print("ERROR: Std. of NVE energy deviation is too large.")
    raise ValueError("Std. of NVE energy deviation is too large.")
    sys.exit(1)
