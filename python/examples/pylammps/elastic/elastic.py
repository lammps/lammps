
from argparse import ArgumentParser
from lammps import PyLammps

def potential(lmp, args):
    """ set up potential and minimization """
    ff_string = ' '
    ff_string = ff_string.join(args.elements) # merge all element string to one string
    lmp.kim("interactions", ff_string)

    # Setup neighbor style
    lmp.neighbor(1.0, "nsq")
    lmp.neigh_modify("once no every 1 delay 0 check yes")

    # Setup minimization style
    lmp.min_style(args.min_style)
    lmp.min_modify("dmax ${dmax} line quadratic")

    # Setup output
    lmp.thermo(1)
    lmp.thermo_style("custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz")
    lmp.thermo_modify("norm no")

    return

def displace(lmp, args, idir):
    """computes the response to a small strain """

    if idir == 1:
        lmp.variable("len0 equal {}".format(lmp.variables["lx0"].value))
    elif idir == 2 or idir == 6:
        lmp.variable("len0 equal {}".format(lmp.variables["ly0"].value))
    else:
        lmp.variable("len0 equal {}".format(lmp.variables["lz0"].value))

    # Reset box and simulation parameters
    lmp.clear()
    lmp.box("tilt large")
    lmp.kim("init", args.kim_model, "metal", "unit_conversion_mode")
    lmp.read_restart("restart.equil")
    lmp.change_box("all triclinic")
    potential(lmp, args)

    # Negative deformation
    lmp.variable("delta equal -${up}*${len0}")
    lmp.variable("deltaxy equal -${up}*xy")
    lmp.variable("deltaxz equal -${up}*xz")
    lmp.variable("deltayz equal -${up}*yz")

    if idir == 1:
        lmp.change_box("all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box")
    elif idir == 2:
        lmp.change_box("all y delta 0 ${delta} yz delta ${deltayz} remap units box")
    elif idir == 3:
        lmp.change_box("all z delta 0 ${delta} remap units box")
    elif idir == 4:
        lmp.change_box("all yz delta ${delta} remap units box")
    elif idir == 5:
        lmp.change_box("all xz delta ${delta} remap units box")
    else:
        lmp.change_box("all xy delta ${delta} remap units box")

    # Relax atoms positions
    lmp.min_style(args.min_style)
    lmp.minimize(args.minimize[0], args.minimize[1], int(args.minimize[2]), int(args.minimize[3]))

    # Obtain new stress tensor
    lmp.variable("pxx1 equal {}".format(lmp.eval("pxx")))
    lmp.variable("pyy1 equal {}".format(lmp.eval("pyy")))
    lmp.variable("pzz1 equal {}".format(lmp.eval("pzz")))
    lmp.variable("pxy1 equal {}".format(lmp.eval("pxy")))
    lmp.variable("pxz1 equal {}".format(lmp.eval("pxz")))
    lmp.variable("pyz1 equal {}".format(lmp.eval("pyz")))

    # Compute elastic constant from pressure tensor
    c1neg = lmp.variables["d1"].value
    c2neg = lmp.variables["d2"].value
    c3neg = lmp.variables["d3"].value
    c4neg = lmp.variables["d4"].value
    c5neg = lmp.variables["d5"].value
    c6neg = lmp.variables["d6"].value

    # Reset box and simulation parameters
    lmp.clear()
    lmp.box("tilt large")
    lmp.kim("init", args.kim_model, "metal", "unit_conversion_mode")
    lmp.read_restart("restart.equil")
    lmp.change_box("all triclinic")
    potential(lmp, args)

    # Positive deformation
    lmp.variable("delta equal ${up}*${len0}")
    lmp.variable("deltaxy equal ${up}*xy")
    lmp.variable("deltaxz equal ${up}*xz")
    lmp.variable("deltayz equal ${up}*yz")

    if idir == 1:
        lmp.change_box("all x delta 0 ${delta} xy delta ${deltaxy} xz delta ${deltaxz} remap units box")
    elif idir == 2:
        lmp.change_box("all y delta 0 ${delta} yz delta ${deltayz} remap units box")
    elif idir == 3:
        lmp.change_box("all z delta 0 ${delta} remap units box")
    elif idir == 4:
        lmp.change_box("all yz delta ${delta} remap units box")
    elif idir == 5:
        lmp.change_box("all xz delta ${delta} remap units box")
    else:
        lmp.change_box("all xy delta ${delta} remap units box")

    # Relax atoms positions
    lmp.min_style(args.min_style)
    lmp.minimize(args.minimize[0], args.minimize[1], int(args.minimize[2]), int(args.minimize[3]))

    # Obtain new stress tensor
    lmp.variable("pxx1 equal {}".format(lmp.eval("pxx")))
    lmp.variable("pyy1 equal {}".format(lmp.eval("pyy")))
    lmp.variable("pzz1 equal {}".format(lmp.eval("pzz")))
    lmp.variable("pxy1 equal {}".format(lmp.eval("pxy")))
    lmp.variable("pxz1 equal {}".format(lmp.eval("pxz")))
    lmp.variable("pyz1 equal {}".format(lmp.eval("pyz")))

    # Compute elasic constant from pressure tensor
    c1pos = lmp.variables["d1"].value
    c2pos = lmp.variables["d2"].value
    c3pos = lmp.variables["d3"].value
    c4pos = lmp.variables["d4"].value
    c5pos = lmp.variables["d5"].value
    c6pos = lmp.variables["d6"].value

    # Combine positive and negative
    lmp.variable("C1{} equal {}".format(idir, 0.5*(c1neg+c1pos)))
    lmp.variable("C2{} equal {}".format(idir, 0.5*(c2neg+c2pos)))
    lmp.variable("C3{} equal {}".format(idir, 0.5*(c3neg+c3pos)))
    lmp.variable("C4{} equal {}".format(idir, 0.5*(c4neg+c4pos)))
    lmp.variable("C5{} equal {}".format(idir, 0.5*(c5neg+c5pos)))
    lmp.variable("C6{} equal {}".format(idir, 0.5*(c6neg+c6pos)))

    return

def elastic():
    """ Compute elastic constant tensor for a crystal

     In order to calculate the elastic constants correctly, care must be taken to specify
     the correct units (units). It is also  important to verify that the minimization of energy
     w.r.t atom  positions in the deformed cell is fully converged.
     One indication of this is that the elastic constants are insensitive
     to the choice of the variable ${up}. Another is to check
     the final max and two-norm forces reported in the log file. If you know
     that minimization is not required, you can set maxiter = 0.0 """

    parser = ArgumentParser(description='A python script to compute elastic properties of bulk materials')

    parser.add_argument("input_data_file", help="The full path & name of the lammps data file.")
    parser.add_argument("kim_model", help="the KIM ID of the interatomic model archived in OpenKIM")
    parser.add_argument("elements", nargs='+', default=['Au'], help="a list of N chemical species, which defines a mapping between atom types in LAMMPS to the available species in the OpenKIM model")
    parser.add_argument("--min_style", default="cg", help="which algorithm will be used for minimization from lammps")
    parser.add_argument("--minimize", type=float, nargs=4, default=[1.0e-4, 1.0e-6, 100, 1000], help="minimization parameters")
    parser.add_argument("--up", type=float, default=1.0e-6, help="the deformation magnitude (in strain units)")
    args = parser.parse_args()

    L = PyLammps()

    L.units("metal")

    # Define the finite deformation size.
    #Try several values to verify that results do not depend on it.
    L.variable("up equal {}".format(args.up))

    # Define the amount of random jiggle for atoms. It prevents atoms from staying on saddle points
    atomjiggle = 1.0e-5

    # metal units, elastic constants in GPa
    cfac = 1.0e-4

    # Define minimization parameters
    L.variable("dmax equal 1.0e-2")

    L.boundary("p", "p", "p") # periodic boundary conditions in all three directions
    L.box("tilt large") # to avoid termination if the final simulation box has a high tilt factor

    # use the OpenKIM model to set the energy interactions
    L.kim("init", args.kim_model, "metal", "unit_conversion_mode")

    L.read_data(args.input_data_file)

    potential(L, args)

    # Need to set mass to something, just to satisfy LAMMPS
    mass_dictionary = {'H': 1.00797, 'He': 4.00260, 'Li': 6.941, 'Be': 9.01218, 'B': 10.81, 'C': 12.011, 'N': 14.0067, 'O': 15.9994, 'F': 18.998403, 'Ne': 20.179, 'Na': 22.98977, 'Mg': 24.305, 'Al': 26.98154, 'Si': 28.0855, 'P': 30.97376, 'S': 32.06, 'Cl': 35.453, 'K': 39.0983, 'Ar': 39.948, 'Ca': 40.08, 'Sc': 44.9559, 'Ti': 47.90, 'V': 50.9415, 'Cr': 51.996, 'Mn': 54.9380, 'Fe': 55.847, 'Ni': 58.70, 'Co': 58.9332, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.72, 'Ge': 72.59, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.80, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.9059, 'Zr': 91.22, 'Nb': 92.9064, 'Mo': 95.94, 'Tc': 98, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.4, 'Ag': 107.868, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.69, 'Sb': 121.75, 'I': 126.9045, 'Te': 127.60, 'Xe': 131.30, 'Cs': 132.9054, 'Ba': 137.33, 'La': 138.9055, 'Ce': 140.12, 'Pr': 140.9077, 'Nd': 144.24, 'Pm': 145, 'Sm': 150.4, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.9254, 'Dy': 162.50, 'Ho': 164.9304, 'Er': 167.26, 'Tm': 168.9342, 'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.85, 'Re': 186.207, 'Os': 190.2, 'Ir': 192.22, 'Pt': 195.09, 'Au': 196.9665, 'Hg': 200.59, 'Tl': 204.37, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 209, 'At': 210, 'Rn': 222, 'Fr': 223, 'Ra': 226.0254, 'Ac': 227.0278, 'Pa': 231.0359, 'Th': 232.0381, 'Np': 237.0482, 'U': 238.029}
    for itype in range(1, len(args.elements)+1):
        L.mass(itype, mass_dictionary.get(args.elements[itype-1], 1.0e-20))

    # Compute initial state at zero pressure
    L.fix(3, "all", "box/relax", "aniso", 0.0)
    L.min_style(args.min_style)
    L.minimize(args.minimize[0], args.minimize[1], int(args.minimize[2]), int(args.minimize[3]))

    L.variable("lx0 equal {}".format(L.eval("lx")))
    L.variable("ly0 equal {}".format(L.eval("ly")))
    L.variable("lz0 equal {}".format(L.eval("lz")))

    # These formulas define the derivatives w.r.t. strain components
    L.variable("d1 equal -(v_pxx1-{})/(v_delta/v_len0)*{}".format(L.eval("pxx"), cfac))
    L.variable("d2 equal -(v_pyy1-{})/(v_delta/v_len0)*{}".format(L.eval("pyy"), cfac))
    L.variable("d3 equal -(v_pzz1-{})/(v_delta/v_len0)*{}".format(L.eval("pzz"), cfac))
    L.variable("d4 equal -(v_pyz1-{})/(v_delta/v_len0)*{}".format(L.eval("pyz"), cfac))
    L.variable("d5 equal -(v_pxz1-{})/(v_delta/v_len0)*{}".format(L.eval("pxz"), cfac))
    L.variable("d6 equal -(v_pxy1-{})/(v_delta/v_len0)*{}".format(L.eval("pxy"), cfac))

    L.displace_atoms("all", "random", atomjiggle, atomjiggle, atomjiggle, 87287, "units box")

    # Write restart
    L.unfix(3)
    L.write_restart("restart.equil")

    for idir in range(1, 7):
        displace(L, args, idir)

    postprocess_and_output(L)
    return

def postprocess_and_output(lmp):
    """Compute the moduli and print everything to screen """

    # Output final values
    c11all = lmp.variables["C11"].value
    c22all = lmp.variables["C22"].value
    c33all = lmp.variables["C33"].value

    c12all = 0.5*(lmp.variables["C12"].value + lmp.variables["C21"].value)
    c13all = 0.5*(lmp.variables["C13"].value + lmp.variables["C31"].value)
    c23all = 0.5*(lmp.variables["C23"].value + lmp.variables["C32"].value)

    c44all = lmp.variables["C44"].value
    c55all = lmp.variables["C55"].value
    c66all = lmp.variables["C66"].value

    c14all = 0.5*(lmp.variables["C14"].value + lmp.variables["C41"].value)
    c15all = 0.5*(lmp.variables["C15"].value + lmp.variables["C51"].value)
    c16all = 0.5*(lmp.variables["C16"].value + lmp.variables["C61"].value)

    c24all = 0.5*(lmp.variables["C24"].value + lmp.variables["C42"].value)
    c25all = 0.5*(lmp.variables["C25"].value + lmp.variables["C52"].value)
    c26all = 0.5*(lmp.variables["C26"].value + lmp.variables["C62"].value)

    c34all = 0.5*(lmp.variables["C34"].value + lmp.variables["C43"].value)
    c35all = 0.5*(lmp.variables["C35"].value + lmp.variables["C53"].value)
    c36all = 0.5*(lmp.variables["C36"].value + lmp.variables["C63"].value)

    c45all = 0.5*(lmp.variables["C45"].value + lmp.variables["C54"].value)
    c46all = 0.5*(lmp.variables["C46"].value + lmp.variables["C64"].value)
    c56all = 0.5*(lmp.variables["C56"].value + lmp.variables["C65"].value)

    # Average moduli for cubic crystals
    c11cubic = (c11all + c22all + c33all)/3.0
    c12cubic = (c12all + c13all + c23all)/3.0
    c44cubic = (c44all + c55all + c66all)/3.0

    bulkmodulus = (c11cubic + 2*c12cubic)/3.0
    shearmodulus1 = c44cubic
    shearmodulus2 = (c11cubic - c12cubic)/2.0
    poisson_ratio = 1.0/(1.0 + c11cubic/c12cubic)

    # print results to screen
    print("=========================================")
    print("Components of the Elastic Constant Tensor")
    print("=========================================")

    print("Elastic Constant C11all = {} GPa".format(c11all))
    print("Elastic Constant C22all = {} GPa".format(c22all))
    print("Elastic Constant C33all = {} GPa".format(c33all))

    print("Elastic Constant C12all = {} GPa".format(c12all))
    print("Elastic Constant C13all = {} GPa".format(c13all))
    print("Elastic Constant C23all = {} GPa".format(c23all))

    print("Elastic Constant C44all = {} GPa".format(c44all))
    print("Elastic Constant C55all = {} GPa".format(c55all))
    print("Elastic Constant C66all = {} GPa".format(c66all))

    print("Elastic Constant C14all = {} GPa".format(c14all))
    print("Elastic Constant C15all = {} GPa".format(c15all))
    print("Elastic Constant C16all = {} GPa".format(c16all))

    print("Elastic Constant C24all = {} GPa".format(c24all))
    print("Elastic Constant C25all = {} GPa".format(c25all))
    print("Elastic Constant C26all = {} GPa".format(c26all))

    print("Elastic Constant C34all = {} GPa".format(c34all))
    print("Elastic Constant C35all = {} GPa".format(c35all))
    print("Elastic Constant C36all = {} GPa".format(c36all))

    print("Elastic Constant C45all = {} GPa".format(c45all))
    print("Elastic Constant C46all = {} GPa".format(c46all))
    print("Elastic Constant C56all = {} GPa".format(c56all))

    print("=========================================")
    print("Average properties for a cubic crystal")
    print("=========================================")

    print("Bulk Modulus = {} GPa".format(bulkmodulus))
    print("Shear Modulus 1 = {} GPa".format(shearmodulus1))
    print("Shear Modulus 2 = {} GPa".format(shearmodulus2))
    print("Poisson Ratio = {}".format(poisson_ratio))

    return

if __name__ == "__main__":
    elastic()
