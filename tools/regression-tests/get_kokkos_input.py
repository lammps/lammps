# This script looks for input scripts under examples/ that have pair/fix/compute styles with KOKKOS support
# and print out separate sets of input lists into 4 files:
#   input-list-pair-kk.txt
#   input-list-fix-kk.txt
#   input-list-compute-kk.txt
#   input-list-misc-kk.txt
# These 4 files will be read in by the regression tester run_tests.py

from argparse import ArgumentParser
import subprocess
import sys

# in_style = fix, pair, compute, angle, bond, min
def get_list(in_style, example_toplevel):
    
    with open(f"input-list-{in_style}-kk.txt", "w") as f:
        # find all the pair styles with the kokkos suffix 
        cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep -v npair | grep {in_style} | grep .cpp"
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
            cmd_str = f"grep -rl '{in_style}_style.*{style}' {example_toplevel}/*/in.* "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list = ' '.join(input_list).split()
            #print(f"There are {len(input_list)} input files that contains {in_style} {style}")
            for input in input_list:
                if input != "":
                    f.write(f"{input}\n")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--examples-top-level", dest="example_toplevel", default="", help="Examples top-level")

    args = parser.parse_args()
    example_toplevel = args.example_toplevel

    with open("input-list-pair-kk.txt", "w") as f:
        # find all the pair styles with the kokkos suffix 
        cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep -v npair | grep pair | grep .cpp"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        kokkos_styles = p.stdout.split('\n')
        style_names = []
        for style in kokkos_styles:
            if style != "":
                # replace "pair_[name]_kokkos.cpp" into "[name]"
                style = style.replace("pair_","")
                style = style.replace("_kokkos.cpp","")
                style = style.replace("_","/")
                style_names.append(style)

        for style in style_names:
            cmd_str = f"grep -rl 'pair_style.*{style}' {example_toplevel}/*/in.* "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list = ' '.join(input_list).split()
            #print(f"There are {len(input_list)} input files that contains pair {style}")
            for input in input_list:
                if input != "":
                    f.write(f"{input}\n")

    with open("input-list-fix-kk.txt", "w") as f:

        # find all the fix styles with the kokkos suffix 
        cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep fix | grep .cpp"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        kokkos_styles = p.stdout.split('\n')
        style_names = []
        for style in kokkos_styles:
            if style != "":
                # replace "fix_[name]_kokkos.cpp" into "[name]"
                style = style.replace("fix_","")
                style = style.replace("_kokkos.cpp","")
                style = style.replace("_","/")
                style_names.append(style)

        for style in style_names:
            cmd_str = f"grep -rl 'fix.*{style}' {example_toplevel}/*/in.* "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list = ' '.join(input_list).split()
            #print(f"There are {len(input_list)} input files that contains fix {style}")
            for input in input_list:
                if input != "":
                    f.write(f"{input}\n")

    with open("input-list-compute-kk.txt", "w") as f:
        # find all the compute styles with the kokkos suffix 
        cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep compute | grep .cpp"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        kokkos_styles = p.stdout.split('\n')
        style_names = []
        for style in kokkos_styles:
            if style != "":
                # replace "compute_[name]_kokkos.cpp" into "[name]"
                style = style.replace("compute_","")
                style = style.replace("_kokkos.cpp","")
                style = style.replace("_","/")
                style_names.append(style)

        for style in style_names:
            cmd_str = f"grep -rl 'compute.*{style}' {example_toplevel}/*/in.* "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list = ' '.join(input_list).split()
            #print(f"There are {len(input_list)} input files that contains compute {style}")
            for input in input_list:
                if input != "":
                    f.write(f"{input}\n")


    with open("input-list-misc-kk.txt", "w") as f:
        
        # find all the angle styles with the kokkos suffix 
        cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep angle | grep .cpp"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        kokkos_styles = p.stdout.split('\n')
        style_names = []
        for style in kokkos_styles:
            if style != "":
                # replace "compute_[name]_kokkos.cpp" into "[name]"
                style = style.replace("angle_","")
                style = style.replace("_kokkos.cpp","")
                style = style.replace("_","/")
                style_names.append(style)

        for style in style_names:
            cmd_str = f"grep -rl 'angle_style.*{style}' {example_toplevel}/*/in.* "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list = ' '.join(input_list).split()
            #print(f"There are {len(input_list)} input files that contains angle {style}")
            for input in input_list:
                if input != "":
                    f.write(f"{input}\n")


        # find all the bond styles with the kokkos suffix 
        cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep bond | grep .cpp"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        kokkos_styles = p.stdout.split('\n')
        style_names = []
        for style in kokkos_styles:
            if style != "":
                # replace "compute_[name]_kokkos.cpp" into "[name]"
                style = style.replace("bond_","")
                style = style.replace("_kokkos.cpp","")
                style = style.replace("_","/")
                style_names.append(style)

        for style in style_names:
            cmd_str = f"grep -rl 'bond_style.*{style}' {example_toplevel}/*/in.* "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list = ' '.join(input_list).split()
            #print(f"There are {len(input_list)} input files that contains bond {style}")
            for input in input_list:
                if input != "":
                    f.write(f"{input}\n")

        # find all the min styles with the kokkos suffix 
        cmd_str = f"ls {example_toplevel}/../src/KOKKOS | grep min | grep .cpp"
        p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
        kokkos_styles = p.stdout.split('\n')
        style_names = []
        for style in kokkos_styles:
            if style != "":
                # replace "compute_[name]_kokkos.cpp" into "[name]"
                style = style.replace("min_","")
                style = style.replace("_kokkos.cpp","")
                style = style.replace("_","/")
                style_names.append(style)

        for style in style_names:
            cmd_str = f"grep -rl 'min_style.*{style}' {example_toplevel}/*/in.* "
            p = subprocess.run(cmd_str, shell=True, text=True, capture_output=True)
            input_list = p.stdout.split('\n')
            input_list = ' '.join(input_list).split()
            #print(f"There are {len(input_list)} input files that contains min {style}")
            for input in input_list:
                if input != "":
                    f.write(f"{input}\n")
