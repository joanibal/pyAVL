"""
    based on autoEditReverse.py and autoEditForward.py in ADflow by G. Kenway
"""


# Import modules
import os
import sys
import re
import argparse


# get input and output directories
parser = argparse.ArgumentParser()
parser.add_argument(
    "input_dir",
    help="Directory of input source files",
    type=str,
)
parser.add_argument(
    "output_dir",
    help="Directory of output source files",
    type=str,
)
parser.add_argument(
    "mode",
    help="Mode of ADflow (reverse or forward)",
    type=str,
    choices=["reverse", "forward"],
)
args = parser.parse_args()

# fmt: off
hand_made_include = "      INCLUDE 'AVL.INC'\n"\
                    "      INCLUDE 'AVL_ad_seeds.inc'\n"
# fmt: on

if args.mode == "reverse":
    # Specify file extension
    file_ext = "_b.f"
    tapenade_include = "INCLUDE 'AVL_b.inc'"
elif args.mode == "forward":
    # Specify file extension
    tapenade_include = "INCLUDE 'AVL_d.inc'"
    file_ext = "_d.f"


print("Directory of input source files  :", args.input_dir)
print("Directory of output source files :", args.output_dir)

for f in os.listdir(args.input_dir):
    if f.endswith(file_ext):
        with open(os.path.join(args.input_dir, f), "r") as fid_org, open(
            os.path.join(args.output_dir, f), "w"
        ) as fid_mod:
            print("\nParsing input file", fid_org.name)

            # read to whole file to string and reposition the pointer
            # at the first byte for future reading
            # all_src = fid_org.read()
            # fid_org.seek(0)

            for line in fid_org:

                if tapenade_include in line:
                    # Insert the hand-written include file
                    fid_mod.write(hand_made_include)
                elif "IMPLICIT NONE" in line:
                    pass
                else:
                    fid_mod.write(line)

        print(" Modified file saved", fid_mod.name)
