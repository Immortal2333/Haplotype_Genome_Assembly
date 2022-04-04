# Genome_assembly
Genome haplotype assembly process and evaluation method. \
Author：Xu Wang niue0501011@gmail.com Data:04/04/2022 \

![process](https://github.com/Immortal2333/Genome_assembly/blob/main/process.jpg?raw=true)

## Data preparation
Pacbio Hifi data: *.ccs.bam, *.ccs.bam.pbi \
Illumina HiC data: *.R1.fastq.gz, *.R2.fastq.gz

## Step by Step
Contigs: \
Step1: Extraction hifi.fastq.zg data \
Step2: HiC data sorting \
Step3: hifiasm haplotype typing

Scaffold: \
Step4:

### Step1：Extraction hifi.fastq.zg data
```smrtlink-tools download
conda install -c hcc smrtlink-tools
```
or
```
wget https://www.pacb.com/support/software-downloads/
```
Using its function:
```
bam2fastq -o <output_name> 1.ccs.bam 2.ccs.bam
```
### Step2：HiC data sorting
```
cat 1.R1.fastq.gz 2.R1.fastq.gz ... > R1.fastq.gz
cat 1.R2.fastq.gz 2.R2.fastq.gz ... > R2.fastq.gz
```
### Step3: hifiasm haplotype typing
hifiasm offical website and usage: https://hifiasm.readthedocs.io/en/latest/hic-assembly.html \
Githup download ZIP: https://github.com/chhylp123/hifiasm , suggested latest version 0.16.1-r375.
```
hifiasm -o <output_name> -t 20 --h1 R1.fastq.gz --h2 R2.fastq.gz HiFi-reads.fq.gz
```
gfa2fa: extraction sequence
```
awk '/^S/{print ">"$2;print $3}' <your.hap1.p_ctg.gfa> > hap1.fa
awk '/^S/{print ">"$2;print $3}' <your.hap2.p_ctg.gfa> > hap2.fa
```


