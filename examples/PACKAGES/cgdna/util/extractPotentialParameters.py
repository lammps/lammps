import sys
import datetime

# Specify path and parameter file name
file_path = '../oxDNA3_sequence_dependent_parameters.txt'

# Specify units: lj or real
units = 'real'

class ConversionFactors:
    def __init__(self, real_units: bool = False):
        if real_units:
            self.unit = 'real'
            self.temperature = 3000.0
            self.energy = 5.961612775922496
            self.length = 8.518
        else:
            self.unit = 'lj'
            self.temperature = 1.0
            self.energy = 1.0
            self.length = 1.0

if units=='lj':
    conversion_factors = ConversionFactors(real_units=False)
if units=='real':
    conversion_factors = ConversionFactors(real_units=True)

# Scan parameter file for strings

lines_fene_delta = []
start_string = 'FENE_DELTA'
exclude_string = 'FENE_DELTA2'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string) and  not stripped_line.startswith(exclude_string):
            part = stripped_line.split(' ')
            lines_fene_delta.append(part[2])
file.close()

lines_fene_r0 = []
start_string = 'FENE_R0'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_fene_r0.append(part[2])
file.close()

lines_excv_sig = []
start_string = 'EXCL_S5'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_excv_sig.append(part[2])
file.close()

lines_excv_rstar = []
start_string = 'EXCL_R5'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_excv_rstar.append(part[2])
file.close()

lines_stk_r0 = []
start_string = 'STCK_R0'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_r0.append(part[2])
file.close()

lines_stk_rc = []
start_string = 'STCK_RC'
exclude_string1 = 'STCK_RCLOW'
exclude_string2 = 'STCK_RCHIGH'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string) and  not stripped_line.startswith(exclude_string1) and  not stripped_line.startswith(exclude_string2):
            part = stripped_line.split(' ')
            lines_stk_rc.append(part[2])
file.close()

lines_stk_rlo = []
start_string = 'STCK_RLOW'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_rlo.append(part[2])
file.close()

lines_stk_rhi = []
start_string = 'STCK_RHIGH'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_rhi.append(part[2])
file.close()

lines_stk_t4a = []
start_string = 'STCK_THETA4_A'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_t4a.append(part[2])
file.close()

lines_stk_t4ts = []
start_string = 'STCK_THETA4_TS'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_t4ts.append(part[2])
file.close()

lines_xstk_r033 = []
start_string = 'CRST_R0_33'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_r033.append(part[2])
file.close()

lines_xstk_rc33 = []
start_string = 'CRST_RC_33'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_rc33.append(part[2])
file.close()

lines_xstk_rlo33 = []
start_string = 'CRST_RLOW_33'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_rlo33.append(part[2])
file.close()

lines_xstk_rhi33 = []
start_string = 'CRST_RHIGH_33'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_rhi33.append(part[2])
file.close()

lines_xstk_r055 = []
start_string = 'CRST_R0_55'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_r055.append(part[2])
file.close()

lines_xstk_rc55 = []
start_string = 'CRST_RC_55'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_rc55.append(part[2])
file.close()

lines_xstk_rlo55 = []
start_string = 'CRST_RLOW_55'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_rlo55.append(part[2])
file.close()

lines_xstk_rhi55 = []
start_string = 'CRST_RHIGH_55'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_rhi55.append(part[2])
file.close()

lines_xstk_t4a33 = []
start_string = 'CRST_THETA4_A_33'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_t4a33.append(part[2])
file.close()

lines_xstk_t4t033 = []
start_string = 'CRST_THETA4_T0_33'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_t4t033.append(part[2])
file.close()

lines_xstk_t4ts33 = []
start_string = 'CRST_THETA4_TS_33'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_t4ts33.append(part[2])
file.close()

lines_xstk_t4a55 = []
start_string = 'CRST_THETA4_A_55'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_t4a55.append(part[2])
file.close()

lines_xstk_t4t055 = []
start_string = 'CRST_THETA4_T0_55'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_t4t055.append(part[2])
file.close()

lines_xstk_t4ts55 = []
start_string = 'CRST_THETA4_TS_55'

# Open the file and read line by line
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_xstk_t4ts55.append(part[2])
file.close()

# Print (converted) potential parameters

date = datetime.datetime.now()

print(("# DATE: %s-%s-%s UNITS: %s CONTRIBUTOR: Oliver Henrich, oliver.henrich@strath.ac.uk CITATION: Bonato, et. al., TBA (2025)") %\
         (date.strftime("%Y"), date.strftime("%m"),date.strftime("%d"),conversion_factors.unit))
print("#")

# Write fene
print("*   fene    %21.15le" % (2.0*conversion_factors.energy), end=" ")
for i in range(len(lines_fene_delta)):
    print("%21.15le" % (float(lines_fene_delta[i])*conversion_factors.length), end=" ")
print(end=" ")
for i in range(len(lines_fene_r0)):
    print("%21.15le" % (float(lines_fene_r0[i])*conversion_factors.length), end=" ")
print("")

# Write excv
print("* * excv    %21.15le %21.15le %21.15le %21.15le %21.15le %21.15le %21.15le %21.15le %21.15le" %\
        (2.0*conversion_factors.energy, 0.7*conversion_factors.energy, 0.675*conversion_factors.energy,\
         2.0*conversion_factors.energy, 0.515*conversion_factors.energy, 0.5*conversion_factors.energy,\
         2.0*conversion_factors.energy, 0.33*conversion_factors.energy, 0.32*conversion_factors.energy), end=" ")
for i in range(len(lines_excv_sig)):
    print("%21.15le" % (float(lines_excv_sig[i])*conversion_factors.length), end=" ")
print(end=" ")
for i in range(len(lines_excv_rstar)):
    print("%21.15le" % (float(lines_excv_rstar[i])*conversion_factors.length), end=" ")
print("")

# Write stk
print("* * stk     %21.15le %21.15le %21.15le" % (1.3523*conversion_factors.energy, \
        2.6717*conversion_factors.energy/conversion_factors.temperature, 6.0/conversion_factors.length), end=" ")
for i in range(len(lines_stk_r0)):
    print("%21.15le" % (float(lines_stk_r0[i])*conversion_factors.length), end=" ")
print(end=" ")
for i in range(len(lines_stk_rc)):
    print("%21.15le" % (float(lines_stk_rc[i])*conversion_factors.length), end=" ")
print(end=" ")
for i in range(len(lines_stk_rlo)):
    print("%21.15le" % (float(lines_stk_rlo[i])*conversion_factors.length), end=" ")
print(end=" ")
for i in range(len(lines_stk_rhi)):
    print("%21.15le" % (float(lines_stk_rhi[i])*conversion_factors.length), end=" ")
print(end=" ")
for i in range(len(lines_stk_t4a)):
    print("%21.15le" % float(lines_stk_t4a[i]), end=" ")
print(end=" ")
print("0", end=" ")
print(end=" ")
for i in range(len(lines_stk_t4ts)):
    print("%21.15le" % float(lines_stk_t4ts[i]), end=" ")
print(end=" ")
print("1.0 0.0 0.901249 1.0 0.0 0.901249 2.0 0.65 2.0 0.65", end=" ")
print("")

# Write hbond
print("* * hbond   %21.15le %21.15le %21.15le %21.15le %21.15le %21.15le 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 1.3 3.141592653589793 0.7904477796209515 2.0 3.141592653589793 0.63727937358744 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45" %\
        (0*conversion_factors.energy, 8.0/conversion_factors.length, 0.4*conversion_factors.length,\
         0.75*conversion_factors.length, 0.34*conversion_factors.length, 0.7*conversion_factors.length), end=" ")
print("")
print("1 4 hbond   %21.15le %21.15le %21.15le %21.15le %21.15le %21.15le 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 1.3 3.141592653589793 0.7904477796209515 2.0 3.141592653589793 0.63727937358744 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45" %\
        (1.0678*conversion_factors.energy, 8.0/conversion_factors.length, 0.4*conversion_factors.length,\
         0.75*conversion_factors.length, 0.34*conversion_factors.length, 0.7*conversion_factors.length), end=" ")
print("")
print("2 3 hbond   %21.15le %21.15le %21.15le %21.15le %21.15le %21.15le 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 1.3 3.141592653589793 0.7904477796209515 2.0 3.141592653589793 0.63727937358744 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45" %\
        (1.0678*conversion_factors.energy, 8.0/conversion_factors.length, 0.4*conversion_factors.length,\
         0.75*conversion_factors.length, 0.34*conversion_factors.length, 0.7*conversion_factors.length), end=" ")
print("")

# Write xstk
print("* * xstk    %21.15le" % (76.0*conversion_factors.energy/conversion_factors.length/conversion_factors.length), end=" ")
for i in range(len(lines_xstk_r033)):
    print("%21.15le" % (float(lines_xstk_r033[i])*conversion_factors.length), end=" ")
for i in range(len(lines_xstk_rc33)):
    print("%21.15le" % (float(lines_xstk_rc33[i])*conversion_factors.length), end=" ")
for i in range(len(lines_xstk_rlo33)):
    print("%21.15le" % (float(lines_xstk_rlo33[i])*conversion_factors.length), end=" ")
for i in range(len(lines_xstk_rhi33)):
    print("%21.15le" % (float(lines_xstk_rhi33[i])*conversion_factors.length), end=" ")
for i in range(len(lines_xstk_r055)):
    print("%21.15le" % (float(lines_xstk_r055[i])*conversion_factors.length), end=" ")
for i in range(len(lines_xstk_rc55)):
    print("%21.15le" % (float(lines_xstk_rc55[i])*conversion_factors.length), end=" ")
for i in range(len(lines_xstk_rlo55)):
    print("%21.15le" % (float(lines_xstk_rlo55[i])*conversion_factors.length), end=" ")
for i in range(len(lines_xstk_rhi55)):
    print("%21.15le" % (float(lines_xstk_rhi55[i])*conversion_factors.length), end=" ")
print("2.25 0.791592653589793 0.58 1.7 1.0 0.68 1.7 1.0 0.68", end=" ")
for i in range(len(lines_xstk_t4a33)):
    print("%21.15le" % float(lines_xstk_t4a33[i]), end=" ")
for i in range(len(lines_xstk_t4t033)):
    print("%21.15le" % float(lines_xstk_t4t033[i]), end=" ")
for i in range(len(lines_xstk_t4ts33)):
    print("%21.15le" % float(lines_xstk_t4ts33[i]), end=" ")
for i in range(len(lines_xstk_t4a55)):
    print("%21.15le" % float(lines_xstk_t4a55[i]), end=" ")
for i in range(len(lines_xstk_t4t055)):
    print("%21.15le" % float(lines_xstk_t4t055[i]), end=" ")
for i in range(len(lines_xstk_t4ts55)):
    print("%21.15le" % float(lines_xstk_t4ts55[i]), end=" ")
print("1.7 0.875 2.266592653589793 0.68 1.7 0.875 2.266592653589793 0.68", end=" ")
print("")

# Write coaxstk
print("* * coaxstk %21.15le %21.15le %21.15le %21.15le %21.15le 2.0 2.891592653589793 0.65 1.3 0.0 0.8 0.9 0.0 0.95 0.9 0.0 0.95 40.0 3.116592653589793" %\
        ((71.53230010095517*conversion_factors.energy/conversion_factors.length/conversion_factors.length),\
         (0.370011*conversion_factors.length), (0.570011*conversion_factors.length),\
         (0.190011*conversion_factors.length), (0.550011*conversion_factors.length)), end=" ")
print("")

# Write dh
print("* * dh      0.815", end=" ")
print("")

# relative stacking strengths
lines_stk_eta_seq = []
lines_stk_eta_value = []

start_string = 'STCK_A_A ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('A_A')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_C_A ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('C_A')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_G_A ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('G_A')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_T_A ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('T_A')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_A_C ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('A_C')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_C_C ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('C_C')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_G_C ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('G_C')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_T_C ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('T_C')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_A_G ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('A_G')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_C_G ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('C_G')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_G_G ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('G_G')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_T_G ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('T_G')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_A_T ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('A_T')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_C_T ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('C_T')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_G_T ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('G_T')
            lines_stk_eta_value.append(part[2])
file.close()

start_string = 'STCK_T_T ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_stk_eta_seq.append('T_T')
            lines_stk_eta_value.append(part[2])
file.close()

print(" ")
print("RELATIVE STACKING AND COAXIAL STACKING STRENGHTS -- REMOVE FROM UPPER POTENTIAL FILE")
for i in range(len(lines_stk_eta_value)):
    print(lines_stk_eta_seq[i], end=" ")
    print(float(lines_stk_eta_value[i])*(1.0-0.18+(0.1*9.0*0.18))/(1.3523+2.6717*0.1))


# relative hydrogen bonding strengths

lines_hb_alpha_seq = []
lines_hb_alpha_value = []

start_string = 'HYDR_A_T ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_hb_alpha_seq.append('A_T')
            lines_hb_alpha_value.append(part[2])
file.close()

start_string = 'HYDR_C_G ='
with open(file_path, 'r') as file:
    for line in file:
        # Remove leading/trailing whitespace
        stripped_line = line.strip()
        # Check if the line starts with the specified string
        if stripped_line.startswith(start_string):
            part = stripped_line.split(' ')
            lines_hb_alpha_seq.append('C_G')
            lines_hb_alpha_value.append(part[2])
file.close()

print(" ")
print("RELATIVE HYDROGEN BONDING STRENGHTS -- REMOVE FROM UPPER POTENTIAL FILE")
for i in range(len(lines_hb_alpha_value)):
    print(lines_hb_alpha_seq[i], end=" ")
    print(float(lines_hb_alpha_value[i])/1.0678)

