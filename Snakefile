import os
os.chdir("/gpfs/data01/bennerlab/home/sjroth/Rat-Pilot/data/")

SAMPLES, = glob_wildcards("raw_data/{sample}_R1.fastq.gz")

rule all:
    input:
        expand("tag_directories/{sample}/tagInfo.txt",sample=SAMPLES)

rule trim:
    input:
        "raw_data/{sample}_R1.fastq.gz"
    output:
        "trimmed/{sample}_R1.fastq.gz"
    threads:
        50
    priority:
        3
    shell:
        "cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -j {threads} -m 20 -o {output} {input}"

rule align:
    input:
        "trimmed/{sample}_R1.fastq.gz"
    output:
        "aligned_files/{sample}.bam"
    threads:
        50
    priority:
        2
    shell:
        """
        bowtie2 -p {threads} -x /gpfs/data01/bennerlab/home/sjroth/software/bowtie_indices/rn6 {input} | samtools sort \
        -@ {threads} -o {output} -
        """

rule make_tag_directory:
    input:
        "aligned_files/{sample}.bam"
    output:
        "tag_directories/{sample}/tagInfo.txt"
    priority:
        1
    shell:
        "makeTagDirectory tag_directories/{wildcards.sample} {input}"