# Diploid_Genome_Assembly
Genome haplotype assembly process and evaluation method. \
Author：Xu Wang niue0501011@gmail.com Data:04/04/2022 

![process](https://github.com/Immortal2333/Diploid_Genome_Assembly/blob/main/process.jpg?raw=true)

## Data preparation
Pacbio Hifi data: *.ccs.bam, *.ccs.bam.pbi \
Illumina HiC data: *.R1.fastq.gz, *.R2.fastq.gz

## Step by Step
Contigs: \
Step1: Extraction hifi.fastq.zg data \
Step2: HiC data sorting \
Step3: hifiasm haplotype typing

Scaffold: \
Step4: Juicer\
Step5: RagTag

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
Usage: hifiasm [options] <in_1.fq> <in_2.fq> <...>

hifiasm -o <output_name> -t 20 --h1 R1.fastq.gz --h2 R2.fastq.gz HiFi-reads.fq.gz
```
gfa2fa: sequence extraction
```
awk '/^S/{print ">"$2;print $3}' <your.hap1.p_ctg.gfa> > hap1.fa
awk '/^S/{print ">"$2;print $3}' <your.hap2.p_ctg.gfa> > hap2.fa
```
### Step4: Juicer
Offical website: https://github.com/aidenlab/juicer \
Download ZIP and install.
```
mkdir juicer && cd juicer
mkdir reference fastq
cat hap1.fa hap2.fa > reference/hap0.fa        # hap1 and hap2 merged to hap0
ln -s R1.fastq.gz fastq/.
ln -s R2.fastq.gz fastq/.
```
Create hap0.fa index
```
bwa index /juicer/reference/hap0.fa
```
Generate site positions
```
Usage: /software/juicer/misc/generate_site_positions.py <restriction enzyme> <out_genome> [location]

python /software/juicer/misc/generate_site_positions.py DpnII hap0 /juicer/reference/hap0.fa
```
Extracte ID and positions
```
awk 'BEGIN{OFS="\t"}{print $1, $NF}' hap0_DpnII.txt > hap0_DpnII.chrom.sizes
```
sh Juicer.sh
```
Usage: juicer.sh [-g genomeID] [-d topDir] [-s site] [-a about] 
                 [-S stage] [-p chrom.sizes path] [-y restriction site file]
                 [-z reference genome file] [-D Juicer scripts directory]
                 [-b ligation] [-t threads] [-h] [-f] [-j]

sh /public/home/wangxu02/software/juicer/CPU/juicer.sh -t 20 -s DpnII -g hap0 \
-d /juicer \
-D /software/juicer \
-z /juicer/reference/hap0.fa \
-p /juicer/hap0_DpnII.chrom.sizes \
-y /juicer/hap0_DpnII.txt
```
Got /juicer/aligned/merged_nodups.txt
### Step5: Ragtag




