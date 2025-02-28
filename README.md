# PipelineProject
Welcome to my pipeline project!

Before you use my wrapper script, there are a few dependencies that you must install (install instructions are for linux os's):
    1.) Bowtie2:
        conda install bowtie2 (install with bioconda)

    2.) SPAdes:
        wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
        tar -xzf SPAdes-4.0.0-Linux.tar.gz
        cd SPAdes-4.0.0-Linux/bin/

    3.) Biopython:
        pip install biopython

    4.) Blast+:
        sudo apt install ncbi-blast+

To call this script follow this format:
python wrapper.py *args*

Required arguments:
    -i1_1 --input1_1: Read 1 in forward direction
    -i1_2 --input1_2: Read 1 in reverse direction
    -i2_1 --input2_1: Read 2 in forward direction
    -i2_2 --input2_2: Read 2 in reverse direction

    -g --genome: The genome that is used to map the reads to for filtering (must be downloaded prior to running)

    -cg --comparing_genomes: The name of the subfamily of genomes to blast against (script will install genomes and create blast library based on this)

Optional Arguments:
    -o --outpath: The path to the directory that the outfile will be created in (defaults to the users home directory)

Included in this repo is sample data which can be used to quickly test the functionality of the pipeline. It can be found at /PipelineProject/sample_data

Example call:
python wrapper.py -g ./PipelineProject/sample_data/RefSeqGenome.fasta -i1_1 ./PipelineProject/sample_data/sampledata1_1.fastq -i1_2 ./PipelineProject/sample_data/sampledata1_2.fastq -2_1 ./PipelineProject/sample_data/sampledata2_1.fastq -i2_2 ./PipelineProject/sample_data/sampledata2_2.fastq -cg Betaherpesvirinae