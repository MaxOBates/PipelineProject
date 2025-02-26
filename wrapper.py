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
    required = True)

    parser.add_argument("-i1_1", "--input1_1",
    help = "input genome 1 read 1 fastq",
    required = True)

    parser.add_argument("-i1_2", "--input1_2",
    help = "input genome 1 read 2 fastq",
    required = True)

    parser.add_argument("-i2_1", "--input2_1",
    help = "input genome 2 read 1 fastq",
    required = True)

    parser.add_argument("-i2_2", "--input2_2",
    help = "input genome 2 read 2 fastq",
    required = True)

    parser.add_argument("-cg", "--comparing_genomes",
    help = "The name of the subfamily of viruses to blast against",
    required = True)

    parser.add_argument("-o", "--outpath",
    help = "the output path to follow, defaults to home directory",
    required = False)

    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
genome_ref = arguments.genome
in_genome_1_1 = arguments.input1_1
in_genome_1_2 = arguments.input1_2
in_genome_2_1 = arguments.input2_1
in_genome_2_2 = arguments.input2_2
comparing_genome = arguments.comparing_genomes
if arguments.outpath:
    out_path = arguments.outpath
else:
    out_path = os.environ["HOME"]

os.chdir(out_path)

try:
    os.system("rm -r PipelineProject_Max_Bates")
except OSError:
    pass

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

def longest_contig(fasta):
    with open(fasta, "r") as f:
        contig = list(SeqIO.parse(f,"fasta"))[0]
        SeqIO.write(contig, "contig_for_blast.fasta", "fasta")

sig_contigs = str(count_sig_contigs("./HCMV-assembly/contigs.fasta"))
sig_bps = str(count_sig_bps("./HCMV-assembly/contigs.fasta"))

with open("PipelineProject.log", 'a') as g:
    g.write(f"There are {sig_contigs} contigs > 1000 bp in the assembly.\n")
    g.write(f"There are {sig_bps} bp in the assembly.\n")



os.system(f"datasets download virus genome taxon {comparing_genome} --include genome")
os.system("unzip ncbi_dataset.zip")
os.system("mv ncbi_dataset/data/genomic.fna ./for_database.fasta")
os.system("rm -r ncbi_dataset ncbi_dataset.zip README.md md5sum.txt")

os.system(f"makeblastdb -in for_database.fasta -out {comparing_genome} -title {comparing_genome} -dbtype nucl")

longest_contig("./HCMV-assembly/contigs.fasta")

os.system(f"blastn -query contig_for_blast.fasta -db {comparing_genome} -out query_out.tsv -outfmt \"6 sacc pident length qstart qend sstart send bitscore evalue stitle\" -max_target_seqs 10 -max_hsps 1")

with open("PipelineProject.log",'a') as f:
    f.write("sacc | pident | length | qstart | qend | sstart | send | bitscore | evalue | stitle\n")
    with open("query_out.tsv", "r") as g:
        blast_out = g.read().splitlines()
    if len(blast_out) <= 10:
        for n in blast_out:
            f.write(n + "\n")
    else:
        for n in range(10):
            f.write(blast_out[n] + "\n")

os.system("mkdir blast_in_out")
os.system(f"mv {comparing_genome}.* for_database.fasta query_out.tsv contig_for_blast.fasta blast_in_out")

os.system("mkdir filtered_reads")
os.system("mv filtered*.*.fq filtered_reads")