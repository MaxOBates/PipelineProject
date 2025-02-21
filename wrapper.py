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
    parser.add_argument("-i1_1", "--input1_1",
    help = "input genome 1 read 1 fastq",
    required = False)
    parser.add_argument("-i1_2", "--input1_2",
    help = "input genome 1 read 2 fastq",
    required = False)
    parser.add_argument("-i2_1", "--input2_1",
    help = "input genome 2 read 1 fastq",
    required = False)
    parser.add_argument("-i2_2", "--input2_2",
    help = "input genome 2 read 2 fastq",
    required = False)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
genome_ref = arguments.genome
in_genome_1_1 = arguments.input1_1
in_genome_1_2 = arguments.input1_2
in_genome_2_1 = arguments.input2_1
in_genome_2_2 = arguments.input2_2

home_dir = os.environ["HOME"]

os.chdir(home_dir)
os.system("mkdir -p PipelineProject_Max_Bates")
os.chdir("./PipelineProject_Max_Bates")

# Ensuring log file is empty
f = open("PipelineProject.log", 'w')
f.close()

### TEST CODE FIX
genome_ref = "~/Github/Genomes/GCF_000845245.1_ViralProj14559_genomic.fna"
in_genome_1_1 = "/home/2025/mbates5/Github/Genomes/SRR5660030_1.fastq"
in_genome_1_2 = "/home/2025/mbates5/Github/Genomes/SRR5660030_2.fastq"
in_genome_2_1 = "/home/2025/mbates5/Github/Genomes/SRR5660033_1.fastq"
in_genome_2_2 = "/home/2025/mbates5/Github/Genomes/SRR5660033_2.fastq"

os.system(f"bowtie2-build {genome_ref} Genome_Index")
os.system(f"bowtie2 --quiet -x Genome_Index -1 {in_genome_1_1} -2 {in_genome_1_2} -S mapped1.sam -p 2 --al-conc-gz filtered1.fq.gz")
os.system(f"bowtie2 --quiet -x Genome_Index -1 {in_genome_2_1} -2 {in_genome_2_2} -S mapped2.sam -p 2 --al-conc-gz filtered2.fq.gz")