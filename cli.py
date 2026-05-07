"""CLI helpers and entry point for the grid-generate command."""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

np.set_printoptions(precision=4, suppress=True, linewidth=100)


# ---------------------------------------------------------------------------
# CLI helpers (formerly in util/util.py)
# ---------------------------------------------------------------------------

def printUsage(NO_ARG_OPTION=False):
    print("Usage is: script.py PATH_TO_URDF (-t FIXED_TARGET_NAMES) (-n FILE_NAMESPACE_NAME) (-d) (-f)")
    print("                    where -d indicates full debug mode")
    print("                    where -f indicates floating base")
    if NO_ARG_OPTION:
        print("Alternative usage assuming grid.cuh is already generated: script.py")


def fileExists(FILE_PATH):
    return pathlib.Path(FILE_PATH).is_file()


def validateFile(FILE_PATH, NO_ARG_OPTION=False):
    if not fileExists(FILE_PATH):
        print("[!Error] grid.cuh does not exist")
        printUsage(NO_ARG_OPTION)
        sys.exit(1)


def parseInputs(NO_ARG_OPTION=False):
    parser = argparse.ArgumentParser(
        description="Process a URDF file and Generate Optimized CUDA Kinematics and Dynamics Code."
    )
    parser.add_argument("urdf_path", nargs="?" if NO_ARG_OPTION else None,
                        help="The path to the URDF file")
    parser.add_argument("-t", "--fixed-target-names", default="", type=str,
                        help="Fixed joint kinematic target names")
    parser.add_argument("-n", "--namespace", default="grid", type=str,
                        help="File namespace name")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Enable debug mode")
    parser.add_argument("-f", "--floating-base", default=False, action="store_true",
                        help="Add a floating base")
    args = parser.parse_args()

    if args.urdf_path is None:
        if NO_ARG_OPTION:
            validateFile("grid.cuh", NO_ARG_OPTION)
            print("Using generated grid.cuh")
            return None
        print("[!Error] No URDF filepath specified")
        printUsage(NO_ARG_OPTION)
        sys.exit(1)

    URDF_PATH = args.urdf_path
    validateFile(URDF_PATH, NO_ARG_OPTION)

    DEBUG_MODE = args.debug
    FLOATING_BASE = args.floating_base
    FILE_NAMESPACE_NAME = args.namespace
    FIXED_TARGET_NAMES = args.fixed_target_names
    if FLOATING_BASE:
        DEBUG_MODE = False

    print("Running with: DEBUG_MODE = " + str(DEBUG_MODE))
    print("           FLOATING_BASE = " + str(FLOATING_BASE))
    print("                    URDF = " + URDF_PATH)
    print("      FIXED_TARGET_NAMES = " + FIXED_TARGET_NAMES)
    print("               FILE_NAME = " + FILE_NAMESPACE_NAME)

    return (URDF_PATH, DEBUG_MODE, FILE_NAMESPACE_NAME, FLOATING_BASE, FIXED_TARGET_NAMES)


def validateRobot(robot, NO_ARG_OPTION=False):
    if robot is None:
        print("[!Error] URDF parsing failed. Please make sure you input a valid URDF file.")
        printUsage(NO_ARG_OPTION)
        sys.exit(1)


# ---------------------------------------------------------------------------
# grid-generate entry point
# ---------------------------------------------------------------------------

def main():
    """Entry point for the ``grid-generate`` CLI command."""
    from URDFParser import URDFParser
    from GRiDCodeGenerator import GRiDCodeGenerator

    URDF_PATH, DEBUG_MODE, FILE_NAMESPACE_NAME, FLOATING_BASE, FIXED_TARGET_NAMES = parseInputs()
    parser = URDFParser()
    robot = parser.parse(URDF_PATH, floating_base=FLOATING_BASE)

    validateRobot(robot)

    codegen = GRiDCodeGenerator(robot, DEBUG_MODE, True, FILE_NAMESPACE=FILE_NAMESPACE_NAME)
    include_homogenous_transforms = not FLOATING_BASE
    codegen.gen_all_code(
        include_homogenous_transforms=include_homogenous_transforms,
        fixed_target_name=FIXED_TARGET_NAMES,
    )
    print("New code generated and saved to grid.cuh!")


if __name__ == "__main__":
    main()
