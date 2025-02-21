import os
import sys
import argparse

#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(
    description="ADD TITLE OF SCRIPT HERE (shows on help -h)")
    parser.add_argument("-g", "--genome",
    help = "ref genome to map to",
    required = False)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
'''infile = arguments.input
outfile = arguments.output'''
genome = arguments.genome

os.system("mkdir -p PipelineProject_Max_Bates")
os.chdir("./PipelineProject_Max_Bates")
os.system("touch PipelineProject.log")
# Ensuring log file is empty
f = open("PipelineProject.log", 'w')
f.close()

### TEST CODE FIX
genome = "~/Github/Genomes/GCF_000845245.1_ViralProj14559_genomic.fna"

os.system(f"bowtie2-build {genome} Genome_Index")

