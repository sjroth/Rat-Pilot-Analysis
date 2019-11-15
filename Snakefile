import os
os.chdir("/gpfs/data01/bennerlab/home/sjroth/Rat-Pilot/atac-data")

SAMPLES,PAIR_ID = glob_wildcards("raw_data/{sample}_{pair_id}.fastq.gz")

rule all:
    input:
        expand("tag_directories/{sample}/tagInfo.txt",sample=SAMPLES),
        expand("aligned_files/{sample}.bam.bai",sample=SAMPLES)

rule trim:
    input:
        "raw_data/{sample}_R1.fastq.gz",
        "raw_data/{sample}_R2.fastq.gz"
    output:
        "trimmed/{sample}_R1.fastq.gz",
        "trimmed/{sample}_R2.fastq.gz"
    threads:
        8
    priority:
        3
    shell:
        """
        cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -j {threads} -m 20 \
        -o {output[0]} -p {output[1]} {input}
        """

rule align:
    input:
        "trimmed/{sample}_R1.fastq.gz",
        "trimmed/{sample}_R2.fastq.gz"
    output:
        "aligned_files/{sample}.bam"
    threads:
        50
    priority:
        2
    shell:
        """
        bowtie2 -p {threads} -x /gpfs/data01/bennerlab/home/sjroth/software/bowtie_indices/rn6 -1 {input[0]} \
        -2 {input[1]} | samtools sort -@ {threads} -o {output} -
        """

rule index_bam:
    input:
        "aligned_files/{sample}.bam"
    output:
        "aligned_files/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule make_tag_directory:
    input:
        "aligned_files/{sample}.bam"
    output:
        "tag_directories/{sample}/tagInfo.txt"
    priority:
        1
    shell:
        "makeTagDirectory tag_directories/{wildcards.sample} {input} -single"