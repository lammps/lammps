# This script looks for input scripts under examples/ that have pair/fix/compute styles with KOKKOS support
# and print out separate sets of input lists into 4 files:
#   input-list-pair-kk.txt
#   input-list-fix-kk.txt
#   input-list-compute-kk.txt
#   input-list-misc-kk.txt
# These 4 files will be read in by the regression tester run_tests.py

from argparse import ArgumentParser
import random
import subprocess
import sys

# in_style = fix, pair, compute, angle, bond, dihedral, improper, min
def generate_list(in_style, example_toplevel, filter_out, output_list):
    
    # find all the pair styles with the kokkos suffix 
    cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep {in_style} | grep .cpp"
    p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
    kokkos_styles = p.stdout.split('\n')
    style_names = []
    for style in kokkos_styles:
        if style != "":
            # replace "{in_style}_[name]_kokkos.cpp" into "[name]"
            style = style.replace(f"{in_style}_","")
            style = style.replace("_kokkos.cpp","")
            style = style.replace("_","/")
            style_names.append(style)

    for style in style_names:
        # find in the in. script a line with "pair_style [name]"
        if in_style == "pair" or in_style == "angle" or in_style == "bond" or in_style == "dihedral" or in_style == "improper":
            cmd_str = f"grep -rl '{in_style}_style.*{style}' {example_toplevel}/*/in.* "
        else:
            # find in the in. script a line with "fix ... [name]" (or "compute ... [name]")
            cmd_str = f"grep -rl '{in_style}.*{style}' {example_toplevel}/*/in.* "

        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        input_list = p.stdout.split('\n')
        input_list = ' '.join(input_list).split()
        for input in input_list:
            if input != "":
                skip = False
                for filter in filter_out:
                    if filter in input:
                        skip = True
                        break
                if skip == True:
                    continue
                else:
                    if input not in output_list:
                        output_list.append(input)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--examples-top-level", dest="example_toplevel", default="", help="Examples top-level")
    parser.add_argument("--filter-out", dest="filter_out", default="", help="Filter out input scripts that contain strings")
    parser.add_argument("--batch-size", dest="batch_size", default=50, help="Batch size of scripts per input list")

    args = parser.parse_args()
    example_toplevel = args.example_toplevel
    filter_out = args.filter_out.split(";")
    batch_size = int(args.batch_size)
    
    # print the list of the input scripts that has each feature to a separate file
    features = [ 'pair', 'fix', 'compute' ]
    for feature in features:
        input_list = []
        generate_list(feature, example_toplevel, filter_out, input_list)

        num_batches = int(len(input_list) / batch_size)
        if num_batches < 2:
            with open(f"input-list-{feature}-kk.txt", "w") as f:
                for input in input_list:
                    if input != "":
                        f.write(f"{input}\n")
        else:
            for idx in range(num_batches):
                with open(f"input-list-{feature}-{idx}-kk.txt", "w") as f:
                    if len(input_list) > batch_size:
                        sampled = random.sample(input_list, batch_size)
                    else:
                        sampled = input_list
                    for input in sampled:
                        if input != "":
                            if input in input_list:
                                input_list.remove(input)
                            f.write(f"{input}\n")

    # combine the list of the input scripts that have these feature to a single file input-list-misc-kk.txt
    features = [ 'angle', 'bond', 'dihedral', 'improper', 'min' ]
    input_list = []
    for feature in features:
        generate_list(feature, example_toplevel, filter_out, input_list)

    with open(f"input-list-misc-kk.txt", "w") as f:
        for input in input_list:
            if input != "":
                f.write(f"{input}\n")

