import os
import sys
import argparse

#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(
    description="ADD TITLE OF SCRIPT HERE (shows on help -h)")
    parser.add_argument("-i", "--input",
    help="input file",
    required=True)
    parser.add_argument("-o", "--output",
    help="output file",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

os.system("mkdir -p PipelineProject_Max_Bates")
os.chdir("./PipelineProject_Max_Bates")
os.system("touch PipelineProject.log")