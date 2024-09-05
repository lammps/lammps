"""
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Kierran Falloon (University of Strathclyde, Glasgow)
                         Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

Program:

Usage: 
$$ python lj2real.py <datafile> <inputfile> [-i]
$$ [-i] flag is optional and is used to convert real -> LJ units.

Requirements:
LAMMPS data file and input file. 

This script assumes a input and data file structure similar to those found in examples/PACKAGES/cgdna/examples/.
"""

import datetime
import os
import sys


class Sections:
    """Sections of the data file"""

    def __init__(
        self,
        bounds: bool,
        masses: bool,
        atoms: bool,
        velocities: bool,
        ellipsoids: bool,
    ):
        self.bounds = bounds  # xlo, xhi, ylo, yhi, zlo, zhi
        self.masses = masses  # Masses
        self.atoms = atoms  # Atoms
        self.velocities = velocities  # Velocities
        self.ellipsoids = ellipsoids  # Ellipsoids


# Conversion factors
class ConversionFactors:
    """Conversion factors for LJ to real units"""

    def __init__(self, invert: bool = False):
        self.inverted = False
        self.temp_conv_factor = 3000.038822
        self.energy_conv_factor = 5.961689060210325
        self.kT_conv_factor = 0.001987204155
        self.mass_conv_factor = 100.0277580236
        self.length_conv_factor = 8.518
        self.time_conv_factor = 1706.0
        self.vel_conv_factor = 0.004992966002344666
        self.angular_mom_conv_factor = 4.254188991883894
        self.density_conv_factor = 0.2687551067436886

        self.oxdna_fene_string = "11.92337812042065 2.1295 6.409795"
        self.oxdna_excv_string = "11.92337812042065 5.9626 5.74965 11.92337812042065 4.38677 4.259 11.92337812042065 2.81094 2.72576"
        self.oxdna_stk_string = "0.70439070204273 3.4072 7.6662 2.72576 6.3885 1.3 0.0 0.8 0.9 0.0 0.95 0.9 0.0 0.95 2.0 0.65 2.0 0.65"
        self.oxdna_hbond_string = "0.0 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592654 0.7 4.0 1.570796327 0.45 4.0 1.570796327 0.45"
        self.oxdna_hbond_1_4_2_3_string = "6.42073911784652 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0 0.7 1.5 0.0 0.7 1.5 0 0.7 0.46 3.141592654 0.7 4.0 1.570796327 0.45 4.0 1.570796327 0.45"
        self.oxdna_xstk_string = "3.9029021145006 4.89785 5.74965 4.21641 5.57929 2.25 0.791592654 0.58 1.7 1.0 0.68 1.7 1.0 0.68 1.5 0 0.65 1.7 0.875 0.68 1.7 0.875 0.68"
        self.oxdna_coaxstk_string = "3.77965257404268 3.4072 5.1108 1.87396 4.94044 2.0 2.541592654 0.65 1.3 0 0.8 0.9 0 0.95 0.9 0 0.95 2 -0.65 2 -0.65"

        self.oxdna2_fene_string = "11.92337812042065 2.1295 6.4430152"
        self.oxdna2_excv_string = "11.92337812042065 5.9626 5.74965 11.92337812042065 4.38677 4.259 11.92337812042065 2.81094 2.72576"
        self.oxdna2_stk_string = "0.70439070204273 3.4072 7.6662 2.72576 6.3885 1.3 0.0 0.8 0.9 0.0 0.95 0.9 0.0 0.95 2.0 0.65 2.0 0.65"
        self.oxdna2_hbond_string = "0.0 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592654 0.7 4.0 1.570796327 0.45 4.0 1.570796327 0.45"
        self.oxdna2_hbond_1_4_2_3_string = "6.36589157849259 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0 0.7 1.5 0.0 0.7 1.5 0 0.7 0.46 3.141592654 0.7 4.0 1.570796327 0.45 4.0 1.570796327 0.45"
        self.oxdna2_xstk_string = "3.9029021145006 4.89785 5.74965 4.21641 5.57929 2.25 0.791592654 0.58 1.7 1.0 0.68 1.7 1.0 0.68 1.5 0 0.65 1.7 0.875 0.68 1.7 0.875 0.68"
        self.oxdna2_coaxstk_string = "4.80673207785863 3.4072 5.1108 1.87396 4.94044 2.0 2.891592653589793 0.65 1.3 0 0.8 0.9 0 0.95 0.9 0 0.95 40.0 3.116592653589793"

        self.oxrna2_fene_string = "11.92337812042065 2.1295 6.482800913"
        self.oxrna2_excv_string = "11.92337812042065 5.9626 5.74965 11.92337812042065 4.38677 4.259 11.92337812042065 2.81094 2.72576"
        self.oxrna2_stk_string = "0.70439070204273 3.66274 7.92174 2.9813 6.64404 0.9 0.0 0.95 0.9 0.0 0.95 1.3 0.0 0.8 1.3 0.0 0.8 2.0 0.65 2.0 0.65"
        self.oxrna2_hbond_string = "0.0 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxrna2_hbond_1_4_2_3_3_4_string = "5.18928666388042 0.93918760272364 3.4072 6.3885 2.89612 5.9626 1.5 0.0 0.7 1.5 0.0 0.7 1.5 0.0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxrna2_xstk_string = "4.92690859644113 4.259 5.1108 3.57756 4.94044 2.25 0.505 0.58 1.7 1.266 0.68 1.7 1.266 0.68 1.7 0.309 0.68 1.7 0.309 0.68"
        self.oxrna2_coaxstk_string = "6.57330882442206 4.259 5.1108 3.57756 4.94044 2.0 2.592 0.65 1.3 0.151 0.8 0.9 0.685 0.95 0.9 0.685 0.95 2.0 -0.65 2.0 -0.65"

        if invert:
            self.invert()

    def invert(self):
        """Inverts the conversion factors for real -> LJ"""
        self.inverted = True
        self.temp_conv_factor = 1.0 / self.temp_conv_factor
        self.energy_conv_factor = 1.0 / self.energy_conv_factor
        self.kT_conv_factor = 1.0 / self.kT_conv_factor
        self.mass_conv_factor = 1.0 / self.mass_conv_factor
        self.length_conv_factor = 1.0 / self.length_conv_factor
        self.time_conv_factor = 1.0 / self.time_conv_factor
        self.vel_conv_factor = 1.0 / self.vel_conv_factor
        self.angular_mom_conv_factor = 1.0 / self.angular_mom_conv_factor
        self.density_conv_factor = 1.0 / self.density_conv_factor

        self.oxdna_fene_string = "2.0 0.25 0.7525"
        self.oxdna_excv_string = "2.0 0.7 0.675 2.0 0.515 0.5 2.0 0.33 0.32"
        self.oxdna_stk_string = (
            "6.0 0.4 0.9 0.32 0.75 1.3 0 0.8 0.9 0 0.95 0.9 0 0.95 2.0 0.65 2.0 0.65"
        )
        self.oxdna_hbond_string = "0.0 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxdna_hbond_1_4_2_3_string = "1.077 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxdna_xstk_string = "47.5 0.575 0.675 0.495 0.655 2.25 0.791592653589793 0.58 1.7 1.0 0.68 1.7 1.0 0.68 1.5 0 0.65 1.7 0.875 0.68 1.7 0.875 0.68"
        self.oxdna_coaxstk_string = "46.0 0.4 0.6 0.22 0.58 2.0 2.541592653589793 0.65 1.3 0 0.8 0.9 0 0.95 0.9 0 0.95 2.0 -0.65 2.0 -0.65"

        self.oxdna2_fene_string = "2.0 0.25 0.7564"
        self.oxdna2_excv_string = "2.0 0.7 0.675 2.0 0.515 0.5 2.0 0.33 0.32"
        self.oxdna2_stk_string = (
            "6.0 0.4 0.9 0.32 0.75 1.3 0 0.8 0.9 0 0.95 0.9 0 0.95 2.0 0.65 2.0 0.65"
        )
        self.oxdna2_hbond_string = "0.0 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxdna2_hbond_1_4_2_3_string = "1.0678 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxdna2_xstk_string = "47.5 0.575 0.675 0.495 0.655 2.25 0.791592653589793 0.58 1.7 1.0 0.68 1.7 1.0 0.68 1.5 0 0.65 1.7 0.875 0.68 1.7 0.875 0.68"
        self.oxdna2_coaxstk_string = "58.5 0.4 0.6 0.22 0.58 2.0 2.891592653589793 0.65 1.3 0 0.8 0.9 0 0.95 0.9 0 0.95 40.0 3.116592653589793"

        self.oxrna2_fene_string = "2.0 0.25 0.761070781051"
        self.oxrna2_excv_string = "2.0 0.7 0.675 2.0 0.515 0.5 2.0 0.33 0.32"
        self.oxrna2_stk_string = "6.0 0.43 0.93 0.35 0.78 0.9 0 0.95 0.9 0 0.95 1.3 0 0.8 1.3 0 0.8 2.0 0.65 2.0 0.65"
        self.oxrna2_hbond_string = "0.0 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxrna2_hbond_1_4_2_3_3_4_string = "0.870439 8.0 0.4 0.75 0.34 0.7 1.5 0 0.7 1.5 0 0.7 1.5 0 0.7 0.46 3.141592653589793 0.7 4.0 1.5707963267948966 0.45 4.0 1.5707963267948966 0.45"
        self.oxrna2_xstk_string = "59.9626 0.5 0.6 0.42 0.58 2.25 0.505 0.58 1.7 1.266 0.68 1.7 1.266 0.68 1.7 0.309 0.68 1.7 0.309 0.68"
        self.oxrna2_coaxstk_string = "80 0.5 0.6 0.42 0.58 2.0 2.592 0.65 1.3 0.151 0.8 0.9 0.685 0.95 0.9 0.685 0.95 2.0 -0.65 2.0 -0.65"


def check_datafile_header(line: str, sections: Sections):
    """Checks for headers to modify corresponding data, since datafile is split into headers.
        Modifies the Sections object to keep track of the current section.

    Args:
        line (str): The line to check
        masses_section (bool): If the current section is the masses section
        atoms_section (bool): If the current section is the atoms section
        velocities_section (bool): If the current section is the velocities section
        ellipsoids_section (bool): If the current section is the ellipsoids section
    """

    if any(header in line for header in ["xlo", "xhi", "ylo", "yhi", "zlo", "zhi"]):
        sections.bounds = True
        sections.masses = False
        sections.atoms = False
        sections.velocities = False
        sections.ellipsoids = False
    elif "Masses" in line:
        sections.bounds = False
        sections.masses = True
        sections.atoms = False
        sections.velocities = False
        sections.ellipsoids = False
    elif "Atoms" in line:
        sections.bounds = False
        sections.masses = False
        sections.atoms = True
        sections.velocities = False
        sections.ellipsoids = False
    elif "Velocities" in line:
        sections.bounds = False
        sections.masses = False
        sections.atoms = False
        sections.velocities = True
        sections.ellipsoids = False
    elif "Ellipsoids" in line:
        sections.bounds = False
        sections.masses = False
        sections.atoms = False
        sections.velocities = False
        sections.ellipsoids = True
    elif "Bonds" in line:
        sections.bounds = False
        sections.masses = False
        sections.atoms = False
        sections.velocities = False
        sections.ellipsoids = False


def modify_datafile(datafile_path: str, conversion_factors: ConversionFactors):
    """Modifies the file by header to use real units.

    Args:
        datafile_path (str): The path to the file to modify
    """
    lines_changed = 0
    current_section = Sections(False, False, False, False, False)

    with open(datafile_path, "r", encoding="UTF-8") as file:
        lines = file.readlines()
        if conversion_factors.inverted:
            lines[0] = (
                "LAMMPS data file in LJ units via oxdna lj2real.py, date "
                + str(datetime.date.today())
                + "\n"
            )
        else:
            lines[0] = (
                "LAMMPS data file in real units via oxdna lj2real.py, date "
                + str(datetime.date.today())
                + "\n"
            )

    for i, line in enumerate(lines):
        check_datafile_header(line, current_section)  # check for headers

        elements = line.split()
        if (
            not elements
            or elements[0] == "#"
            or any(
                header in line
                for header in ["Masses", "Atoms", "Velocities", "Ellipsoids", "Bonds"]
            )
        ):
            continue

        # modify the line based on the current section it is in
        if current_section.bounds:
            elements[0:2] = [
                str(int(float(x) * conversion_factors.length_conv_factor))
                for x in elements[0:2]
            ]
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1
        if current_section.masses:
            elements[1] = str(
                round(float(elements[1]) * conversion_factors.mass_conv_factor, 4)
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1
        elif current_section.atoms:
            elements[2:5] = [
                str(float(x) * conversion_factors.length_conv_factor)
                for x in elements[2:5]
            ]
            elements[7] = str(
                float(elements[7]) * conversion_factors.density_conv_factor
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1
        elif current_section.velocities:
            elements[1:4] = [
                str(float(x) * conversion_factors.vel_conv_factor)
                for x in elements[1:4]
            ]
            elements[4:7] = [
                str(float(x) * conversion_factors.angular_mom_conv_factor)
                for x in elements[4:7]
            ]
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1
        elif current_section.ellipsoids:
            elements[1:4] = [
                str(float(x) * conversion_factors.length_conv_factor)
                for x in elements[1:4]
            ]
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

    if conversion_factors.inverted:
        new_datafile_path = datafile_path + "_lj"
    else:
        new_datafile_path = datafile_path + "_real"

    with open(new_datafile_path, "w", encoding="UTF-8") as file:
        file.writelines(lines)
        if lines_changed == 0:
            print(
                "Warning: No lines changed in data file. Ensure correct usage: python lj2real.py <datafile> <inputfile> [-i]"
            )
        else:
            print(f"Data file lines changed: {lines_changed}")

    return new_datafile_path


def modify_inputfile(inputfile_path: str, conversion_factors: ConversionFactors):
    """Modifies the input file line by line to use real units.

    Args:
        inputfile_path (str): The path to the input file to modify
    """

    lines_changed = 0
    oxdna2_flag, oxrna2_flag = False, False

    with open(inputfile_path, "r", encoding="UTF-8") as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if "oxdna2" in line and not oxdna2_flag:
            oxdna2_flag = True
            print("Note: oxdna2 found in input file. Using oxdna2 conversion factors.")
        if "oxrna2" in line and not oxrna2_flag:
            oxrna2_flag = True
            print("Note: oxrna2 found in input file. Using oxrna2 conversion factors.")
        if oxdna2_flag and oxrna2_flag:
            print(
                "Warning: Both oxdna2 and oxrna2 found in input file. Output will likely be incorrect."
            )

        if "variable T" in line:
            old_value = line.split()[3]

            new_value = str(
                round(float(old_value) * conversion_factors.temp_conv_factor, 1)
            )
            lines[i] = line.replace(old_value, new_value)
            lines_changed += 1

        elif "units" in line:
            if conversion_factors.inverted:
                lines[i] = "units lj\n"
            else:
                lines[i] = "units real\n"
            lines_changed += 1

        elif "atom_modify" in line:
            elements = line.split()
            elements[3] = str(
                round(float(elements[3]) * conversion_factors.length_conv_factor, 3)
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

        elif "neighbor" in line:
            elements = line.split()
            elements[1] = str(
                round(float(elements[1]) * conversion_factors.length_conv_factor, 3)
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

        elif "read_data" in line:
            elements = line.split()
            if conversion_factors.inverted:
                elements[1] = elements[1] + "_lj"
            else:
                elements[1] = (
                    elements[1] + "_real"
                )  # naming convention of datafile after conversion
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

        elif "mass" in line:
            elements = line.split()
            elements[4] = str(
                round(float(elements[4]) * conversion_factors.mass_conv_factor, 4)
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

        elif "bond_coeff" in line or "pair_coeff" in line:
            if ".lj" in line or ".real" in line:
                if conversion_factors.inverted:
                    line = line.replace(".real", ".lj")
                else:
                    line = line.replace(".lj", ".real")
                lines[i] = line
                lines_changed += 1

                if "stk" in line and "xstk" not in line and "coaxstk" not in line:
                    elements = line.split()
                    elements[6] = str(  # convert xi
                        round(
                            float(elements[6]) * conversion_factors.energy_conv_factor,
                            14,
                        )
                    )
                    elements[7] = str(  # convert kappa
                        round(float(elements[7]) * conversion_factors.kT_conv_factor, 9)
                    )
                    lines[i] = " ".join(elements) + "\n"

            else:
                elements = line.split()

                if "bond_coeff" in line:
                    if oxdna2_flag:
                        elements[2:] = conversion_factors.oxdna2_fene_string.split()
                    elif oxrna2_flag:
                        elements[2:] = conversion_factors.oxrna2_fene_string.split()
                    else:
                        elements[2:] = conversion_factors.oxdna_fene_string.split()

                elif "excv" in line:
                    if oxdna2_flag:
                        elements[4:] = conversion_factors.oxdna2_excv_string.split()
                    elif oxrna2_flag:
                        elements[4:] = conversion_factors.oxrna2_excv_string.split()
                    else:
                        elements[4:] = conversion_factors.oxdna_excv_string.split()

                elif "stk" in line:

                    if "coaxstk" in line:
                        if oxdna2_flag:
                            elements[4:] = (
                                conversion_factors.oxdna2_coaxstk_string.split()
                            )
                        elif oxrna2_flag:
                            elements[4:] = (
                                conversion_factors.oxrna2_coaxstk_string.split()
                            )
                        else:
                            elements[4:] = (
                                conversion_factors.oxdna_coaxstk_string.split()
                            )

                    elif "xstk" in line:
                        if oxdna2_flag:
                            elements[4:] = conversion_factors.oxdna2_xstk_string.split()
                        elif oxrna2_flag:
                            elements[4:] = conversion_factors.oxrna2_xstk_string.split()
                        else:
                            elements[4:] = conversion_factors.oxdna_xstk_string.split()

                    else:  # stk
                        elements[6] = str(  # convert xi
                            round(
                                float(elements[6])
                                * conversion_factors.energy_conv_factor,
                                14,
                            )
                        )
                        elements[7] = str(  # convert kappa
                            round(
                                float(elements[7]) * conversion_factors.kT_conv_factor,
                                9,
                            )
                        )
                        if oxdna2_flag:
                            elements[8:] = conversion_factors.oxdna2_stk_string.split()
                        elif oxrna2_flag:
                            elements[8:] = conversion_factors.oxrna2_stk_string.split()
                        else:
                            elements[8:] = conversion_factors.oxdna_stk_string.split()

                elif "hbond" in line:
                    if elements[1] == "*" and elements[2] == "*":
                        if oxdna2_flag:
                            elements[5:] = (
                                conversion_factors.oxdna2_hbond_string.split()
                            )
                        elif oxrna2_flag:
                            elements[5:] = (
                                conversion_factors.oxrna2_hbond_string.split()
                            )
                        else:
                            elements[5:] = conversion_factors.oxdna_hbond_string.split()
                    else:
                        if oxdna2_flag:
                            elements[5:] = (
                                conversion_factors.oxdna2_hbond_1_4_2_3_string.split()
                            )
                        elif oxrna2_flag:
                            elements[5:] = (
                                conversion_factors.oxrna2_hbond_1_4_2_3_3_4_string.split()
                            )
                        else:
                            elements[5:] = (
                                conversion_factors.oxdna_hbond_1_4_2_3_string.split()
                            )
                lines[i] = " ".join(elements) + "\n"
                lines_changed += 1

        elif "langevin" in line:
            elements = line.split()
            elements[6] = str(
                round(float(elements[6]) * conversion_factors.time_conv_factor, 2)
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

        elif "timestep" in line:
            elements = line.split()
            elements[1] = str(
                round(float(elements[1]) * conversion_factors.time_conv_factor, 5)
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

        elif "comm_modify" in line:
            elements = line.split()
            elements[2] = str(
                round(float(elements[2]) * conversion_factors.length_conv_factor, 1)
            )
            lines[i] = " ".join(elements) + "\n"
            lines_changed += 1

        else:
            continue

    if conversion_factors.inverted:
        new_inputfile_path = inputfile_path + "_lj"
    else:
        new_inputfile_path = inputfile_path + "_real"

    with open(new_inputfile_path, "w", encoding="UTF-8") as file:
        if conversion_factors.inverted:
            file.write(
                "# LAMMPS input file in LJ units via oxdna lj2real.py, date "
                + str(datetime.date.today())
                + "\n"
            )
        else:
            file.write(
                "# LAMMPS input file in real units via oxdna lj2real.py, date "
                + str(datetime.date.today())
                + "\n"
            )
        file.writelines(lines)
        if lines_changed == 0:
            print(
                "Warning: No lines changed in input file. Ensure correct usage: python lj2real.py <datafile> <inputfile> [-i]"
            )
        else:
            print(f"Input file lines changed: {lines_changed}")

    return new_inputfile_path


def main():
    """Main function"""

    print(
        "\nLAMMPS data and input file conversion to real units via oxdna convert_data.py\n"
        "Note: This script assumes a input and data file structure similar to those found in examples/PACKAGES/cgdna/examples/.\n"
        "Ensure output is checked for correctness."
    )

    if len(sys.argv) > 2:
        datafile_path = sys.argv[1]
        inputfile_path = sys.argv[2]
        invert = False
        if len(sys.argv) > 3:
            if sys.argv[3] == "-i":
                invert = True
                print("Performing real -> LJ conversion.")
            else:
                print(
                    "Invalid flag. Usage: python lj2real.py <datafile> <inputfile> [-i]"
                )
                sys.exit(1)

        if invert:
            conversion_factors = ConversionFactors(invert=True)
            print(
                "\nUsing conversion factors (real T, m, l, t, v, ρ -> LJ T*, m*, l*, t*, v*, ρ*):"
            )

        else:
            conversion_factors = ConversionFactors(invert=False)
            print(
                "\nUsing conversion factors (LJ T*, m*, l*, t*, v*, ρ* -> real T, m, l, t, v, ρ):"
            )

    else:
        print("\nUsage: python lj2real.py <datafile> <inputfile> [-i]")
        print("\t[-i] flag is optional and is used to convert real -> LJ units.\n")
        sys.exit(1)

    conversion_factors_string = (
        f"Temperature T\t≈ {round(conversion_factors.temp_conv_factor, 6)} T* (K)\n"
        f"Energy ε\t≈ {round(conversion_factors.energy_conv_factor, 6)} ε* (kcal/mol)\n"
        f"Mass m\t\t≈ {round(conversion_factors.mass_conv_factor, 6)} m* (g/mol)\n"
        f"Length l\t≈ {round(conversion_factors.length_conv_factor, 6)} l* (Å)\n"
        f"Time t\t\t≈ {round(conversion_factors.time_conv_factor, 6)} t* (fs)\n"
        f"Velocity v\t≈ {round(conversion_factors.vel_conv_factor, 6)} v* (Å/fs)\n"
        f"AngMom l\t≈ {round(conversion_factors.angular_mom_conv_factor, 6)} (g/mol Å^2/fs)\n"
        f"Density ρ\t≈ {round(conversion_factors.density_conv_factor, 6)} ρ* (g/cm^3)\n"
        f"Calculated using Sengar, Ouldridge, Henrich, Rovigatti, & Šulc. Front Mol Biosci 8 (2021). & https://docs.lammps.org/units.html.\n"
        f"See examples/PACKAGES/cgdna/util/lj2real.py for exact conversion factors.\n"
    )
    print(conversion_factors_string)

    print("Current directory: ", os.getcwd())

    try:
        new_datafile_path = modify_datafile(datafile_path, conversion_factors)
        print(f"New data file: {new_datafile_path}")
    except Exception as e:
        print(f"Error modifying data file: {e}")

    try:
        new_inputfile_path = modify_inputfile(inputfile_path, conversion_factors)
        print(f"New input file: {new_inputfile_path}\n")
    except Exception as e:
        print(f"Error modifying input file: {e}")


if __name__ == "__main__":
    main()
