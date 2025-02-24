import os
import sys
import argparse
from Bio import SeqIO

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

def read_count(fastq_file_list):
    read_counter = 0
    for fastq in fastq_file_list:
        with open(fastq) as f:
            for record in SeqIO.parse(f,"fastq"):
                read_counter += 1
    return read_counter

num_reads_before1 = read_count([in_genome_1_1,in_genome_1_2])
num_reads_before2 = read_count([in_genome_2_1,in_genome_2_2])

os.system(f"bowtie2-build {genome_ref} Genome_Index")
os.system(f"bowtie2 --quiet -x Genome_Index -1 {in_genome_1_1} -2 {in_genome_1_2} -S mapped1.sam -p 2 --al-conc filtered1.fq")
os.system(f"bowtie2 --quiet -x Genome_Index -1 {in_genome_2_1} -2 {in_genome_2_2} -S mapped2.sam -p 2 --al-conc filtered2.fq")

os.system("mkdir Genome_Index")
os.system("mv Genome_Index.* Genome_Index")

os.system("mkdir Mapped_Reads")
os.system("mv mapped*.sam Mapped_Reads")

num_reads_after1 = read_count(["filtered1.1.fq","filtered1.2.fq"])
num_reads_after2 = read_count(["filtered2.1.fq","filtered2.2.fq"])

with open("PipelineProject.log", 'a') as f:
    f.write(f"Donor1 (2dpi) had {num_reads_before1} read pairs before Bowtie2 filtering and {num_reads_after1} read pairs after.\n")
    f.write(f"Donor2 (6dpi) had {num_reads_before2} read pairs before Bowtie2 filtering and {num_reads_after2} read pairs after.\n")

spades_cmd = "spades.py -k 99 -t 2 --only-assembler --pe1-1 filtered1.1.fq --pe1-2 filtered1.2.fq --pe2-1 filtered2.1.fq --pe2-2 filtered2.2.fq -o HCMV-assembly/"
os.system(spades_cmd)

with open("PipelineProject.log", 'a') as f:
    f.write(spades_cmd + "\n")

def count_sig_contigs(fasta):
    contig_counter = 0
    with open(fasta, "r") as f:    
        for record in SeqIO.parse(f,"fasta"):
            if len(record.seq) > 1000:
                contig_counter += 1
    return contig_counter

def count_sig_bps(fasta):
    bp_counter = 0
    with open(fasta, "r") as f:    
        for record in SeqIO.parse(f,"fasta"):
            if len(record.seq) > 1000:
                bp_counter += len(record.seq)
    return bp_counter

sig_contigs = str(count_sig_contigs("./HCMV-assembly/contigs.fasta"))
sig_bps = str(count_sig_bps("./HCMV-assembly/contigs.fasta"))

with open("PipelineProject.log", 'a') as g:
    g.write(f"There are {sig_contigs} contigs > 1000 bp in the assembly.\n")
    g.write(f"There are {sig_bps} bp in the assembly.\n")
