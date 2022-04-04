# Genom_assembly
Genome haplotype assembly process and evaluation method.

# Data preparation
Pacbio Hifi data: *.ccs.bam, *.ccs.bam.pbi \
Illumina HiC data: *.R1.fastq.gz, *.R2.fastq.gz

Step1：bam2fastq

conda install -c hcc smrtlink-tools

bam2fastq -o hifiasm.gz m001.bam m002.bam

Step2：HiC data
