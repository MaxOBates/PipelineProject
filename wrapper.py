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

# Ensuring that the out directory does not exist
try:
    os.system("rm -r PipelineProject_Max_Bates")
except OSError:
    pass

# Making out directory and making it the working directory
os.system("mkdir -p PipelineProject_Max_Bates")
os.chdir("./PipelineProject_Max_Bates")

# Ensuring out log file is empty
f = open("PipelineProject.log", 'w')
f.close()

# finding the number of reads in infile
def read_count(fastq_file_list):
    read_counter = 0
    for fastq in fastq_file_list:
        with open(fastq) as f:
            for record in SeqIO.parse(f,"fastq"):
                read_counter += 1
    return read_counter

# calling readcount for all infiles
num_reads_before1 = read_count([in_genome_1_1])
num_reads_before2 = read_count([in_genome_2_1])

# running bowtie2 against genome to compare to
os.system(f"bowtie2-build {genome_ref} Genome_Index")
os.system(f"bowtie2 --quiet -x Genome_Index -1 {in_genome_1_1} -2 {in_genome_1_2} -S mapped1.sam -p 2 --al-conc filtered1.fq")
os.system(f"bowtie2 --quiet -x Genome_Index -1 {in_genome_2_1} -2 {in_genome_2_2} -S mapped2.sam -p 2 --al-conc filtered2.fq")

# cleaning up
os.system("mkdir Genome_Index")
os.system("mv Genome_Index.* Genome_Index")
os.system("mkdir Mapped_Reads")
os.system("mv mapped*.sam Mapped_Reads")

# finding the number of reads after filtering
num_reads_after1 = read_count(["filtered1.1.fq"])
num_reads_after2 = read_count(["filtered2.1.fq"])

# writing number of reads before and after filtering to log file
with open("PipelineProject.log", 'a') as f:
    f.write(f"Donor1 (2dpi) had {num_reads_before1} read pairs before Bowtie2 filtering and {num_reads_after1} read pairs after.\n")
    f.write(f"Donor2 (6dpi) had {num_reads_before2} read pairs before Bowtie2 filtering and {num_reads_after2} read pairs after.\n")

# calling spades
spades_cmd = "spades.py -k 99 -t 2 --only-assembler --pe1-1 filtered1.1.fq --pe1-2 filtered1.2.fq --pe2-1 filtered2.1.fq --pe2-2 filtered2.2.fq -o HCMV-assembly/"
os.system(spades_cmd)

# writing the spades command to log file
with open("PipelineProject.log", 'a') as f:
    f.write(spades_cmd + "\n")

# counting contigs > 1000 bp
def count_sig_contigs(fasta):
    contig_counter = 0
    with open(fasta, "r") as f:    
        for record in SeqIO.parse(f,"fasta"):
            if len(record.seq) > 1000:
                contig_counter += 1
    return contig_counter

# counting number of bps in contigs > 1000 bp
def count_sig_bps(fasta):
    bp_counter = 0
    with open(fasta, "r") as f:    
        for record in SeqIO.parse(f,"fasta"):
            if len(record.seq) > 1000:
                bp_counter += len(record.seq)
    return bp_counter

# writing the longest contig to a fasta file
def longest_contig(fasta):
    with open(fasta, "r") as f:
        contig = list(SeqIO.parse(f,"fasta"))[0]
        SeqIO.write(contig, "contig_for_blast.fasta", "fasta")

# calling count_sig_contigs/_bps
sig_contigs = str(count_sig_contigs("./HCMV-assembly/contigs.fasta"))
sig_bps = str(count_sig_bps("./HCMV-assembly/contigs.fasta"))

# writing output of count_sig_contigs to log file
with open("PipelineProject.log", 'a') as g:
    g.write(f"There are {sig_contigs} contigs > 1000 bp in the assembly.\n")
    g.write(f"There are {sig_bps} bp in the assembly.\n")

# downloading the genome to blast against
os.system(f"datasets download virus genome taxon {comparing_genome} --include genome")
os.system("unzip ncbi_dataset.zip")

# cleaning up
os.system("mv ncbi_dataset/data/genomic.fna ./for_database.fasta")
os.system("rm -r ncbi_dataset ncbi_dataset.zip README.md md5sum.txt")

# making database from download to blast against
os.system(f"makeblastdb -in for_database.fasta -out {comparing_genome} -title {comparing_genome} -dbtype nucl")

# calling longest_contig to blast against local database
longest_contig("./HCMV-assembly/contigs.fasta")
os.system(f"blastn -query contig_for_blast.fasta -db {comparing_genome} -out query_out.tsv -outfmt \"6 sacc pident length qstart qend sstart send bitscore evalue stitle\" -max_target_seqs 10 -max_hsps 1")

# writing top 10 blast results with header to log file
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

# cleaning up
os.system("mkdir blast_in_out")
os.system(f"mv {comparing_genome}.* for_database.fasta query_out.tsv contig_for_blast.fasta blast_in_out")
os.system("mkdir filtered_reads")
os.system("mv filtered*.*.fq filtered_reads")