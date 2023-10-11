# Haplotype Genome Assembly and Evaluation
![Update](https://img.shields.io/badge/Update-10/10/2023-green?logo=github)
![Author](https://img.shields.io/badge/Author-Xu.Wang-orange)
![Email](https://img.shields.io/badge/Email-571720850@qq.com-blue?)

## Data preparation
**Input Data**
* Pacbio HiFi data: `.ccs.bam` or `.fastq.gz`
* Illumina Hi-C data: `.R1.fastq.gz` and `.R2.fastq.gz`
* ONT data: `.bam` or `.fastq.gz` (if you have)

**Dependency** <img src="https://github.com/zhouyflab/Polygenetic_Basis_Seedless_Grapes/blob/main/genome.assembly.jpg" align="right" width="50%">

The detials of all tools can be available in their **offical website** as followed and most of them can quickly install using [Anaconda](https://anaconda.org/):
* [Smrtlink-tools/bam2fastq](https://anaconda.org/hcc/smrtlink-tools) 
* [Fastp](https://github.com/OpenGene/fastp)
* [Hifiasm](https://github.com/chhylp123/hifiasm)
* [Juicer](https://github.com/aidenlab/juicer)
* [RagTag](https://github.com/malonge/RagTag)
* [Juicebox_scripts](https://github.com/phasegenomics/juicebox_scripts)
* [3d-dna](https://github.com/aidenlab/3d-dna)
* [Juicebox](https://github.com/aidenlab/Juicebox/wiki/Download)

## Genome Assembly
### Step0: Data Preprocessing
* If you have raw sequencing bam, please convert `bam` to `fastq` format.
```
bam2fastq -o output.fastq.gz m001.bam m002.bam ...
```
* The illumina Hi-C data need to filter based on quality control (QC).
```
cat file1 file2 .. > R1.fastq.gz
cat file3 file4 .. > R2.fastq.gz
fastp -i R1.fastq.gz -o R1.clean.fastq.gz -I R2.fastq.gz -O R2.clean.fastq.gz -w 6 --html clean.html --json clean.json
```

### Step1：[hifiasm](https://github.com/chhylp123/hifiasm): Hi-C integration assembly
* More detial can be found [hifiasm parameter](https://github.com/chhylp123/hifiasm).
```
hifiasm -o <output_name> -t 20 --h1 R1.fastq.gz --h2 R2.fastq.gz hifiasm.fastq.gz
```
* gfa2fa: sequence extraction
```
awk '/^S/{print ">"$2;print $3}' <your.hap1.p_ctg.gfa> > hap1.fa
awk '/^S/{print ">"$2;print $3}' <your.hap2.p_ctg.gfa> > hap2.fa
```

### Step2: [Juicer](https://github.com/aidenlab/juicer)
* Preparation 
```
mkdir juicer && cd juicer
mkdir reference fastq
cat hap1.fa hap2.fa > reference/hap0.fa        # hap1 and hap2 merged to hap0 in reference file
ln -s R1.clean.fastq.gz fastq/.      # link Hi-C clean data into fastq file
ln -s R2.clean.fastq.gz fastq/.
```
* Create reference index
```
bwa index reference/hap0.fa
```
* Generate site positions. The script, generate_site_positions.py, can be downloaded from **Offical website**.
```
Usage: /your_download_path/juicer/misc/generate_site_positions.py <restriction enzyme> <out_genome> [location]

python /your_download_path/juicer/misc/generate_site_positions.py DpnII hap0 reference/hap0.fa
```
* Extracte ID and positions
```
awk 'BEGIN{OFS="\t"}{print $1, $NF}' hap0_DpnII.txt > hap0_DpnII.chrom.sizes
```
* Run Juicer.sh
```
Usage: juicer.sh [-g genomeID] [-d topDir] [-s site] [-a about] 
                 [-S stage] [-p chrom.sizes path] [-y restriction site file]
                 [-z reference genome file] [-D Juicer scripts directory]
                 [-b ligation] [-t threads] [-h] [-f] [-j]

sh /your_download_path/juicer/CPU/juicer.sh -t 20 -s DpnII -g hap0 \
-d ./ \
-D /software/juicer \
-z reference/hap0.fa \
-p hap0_DpnII.chrom.sizes \
-y hap0_DpnII.txt
```
* Finally, we will get /juicer/aligned/**merged_nodups.txt**

### Step3: [Ragtag](https://github.com/malonge/RagTag/wiki/scaffold)
Scaffolding is the process of ordering and orienting draft assembly (query) sequences into longer sequences. Gaps (stretches of "N" characters) are placed between adjacent query sequences to indicate the presence of unknown sequence. RagTag uses whole-genome alignments to a reference assembly to scaffold query sequences. RagTag does not alter input query sequence in any way and only orders and orients sequences, joining them with gaps.
* Preparation
```
mkdir ragtag && cd ragtag
ln -s hap1.fa ./
ln -s hap2.fa ./
ln -s reference.fa ./       # This reference genome refers to the genome that is recognized as the best assembled in the industry
```
* Run RagTag
```
usage: ragtag scaffold <reference.fa> <query.fa>

python /your_download_path/ragtag/bin/ragtag.py scaffold \
reference.fa hap1.fa -o hap1

python /your_download_path/ragtag/bin/ragtag.py scaffold \
reference.fa hap2.fa -o hap2
```
* agp2assembly.py ([juicebox_script](https://github.com/phasegenomics/juicebox_scripts))
```
usage:  agp2assembly.py <input_agp_file> <output_assembly_file>

python /your_download_path/juicebox_scripts/agp2assembly.py \
hap1/ragtag.scaffold.agp hap1.assembly

python /your_download_path/juicebox_scripts/agp2assembly.py \
hap2/ragtag.scaffold.agp hap2.assembly
```
* Combined the recording of contigs assembly. **The script has uploaded in this pipeline.**
```
python 01_2assemblyto2.py hap1.assembly hap2.assembly hap0.assembly       
```
### Step4: [3d-dna](https://github.com/aidenlab/3d-dna)
```
USAGE: ./run-assembly-visualizer.sh [options] <path_to_input_assembly_file> <path_to_input_mnd_file>

OPTIONS:
    -q mapq                        Build map for a specific mapq threshold (default is 1). 0 is loose, and1 it strict.

/your_download_path/3d-dna/visualize/run-assembly-visualizer.sh -q 0 \
ragtag/hap0.assembly aligned/merged_nodups.txt
```
Notes: Alignment score is the quality of of matching between the read-sequence and reference-sequence. Mapping quality is the confidence that the read is correctly mapped to the genomic coordinates. For example, a read may be mapped to several genomic locations with almost a perfect match in all locations. In that case, alignment score will be high but mapping quality will be low.

### Step5: Manual adjusted genome ([Juicebox](https://github.com/aidenlab/Juicebox/wiki/Download))
In the output of step4 section, we will get two file: `final.assembly` and `0.hic`. Loading these files to **Juicebox** for visualization. 
* [Juicebox Official Website User Manual](https://github.com/aidenlab/Juicebox/wiki)
* [Tutorial video for Juicebox Assembly Tools](https://www.youtube.com/watch?v=Nj7RhQZHM18)
<div style="display:flex;">
  <img1 src="https://github.com/Immortal2333/Haplotype_Genome_Assembly/blob/main/Juicebox_pics/hic.heatmap.png" style="flex:50%; padding:5px;">
  <img2 src="https://github.com/Immortal2333/Haplotype_Genome_Assembly/blob/main/Juicebox_pics/contigs.details.png" style="flex:50%; padding:5px;">
</div>
























