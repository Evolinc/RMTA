[![Docker Pulls](https://img.shields.io/docker/pulls/evolinc/rmta.svg)](https://hub.docker.com/r/evolinc/rmta/)
[![Docker Stars](https://img.shields.io/docker/stars/evolinc/rmta.svg)](https://hub.docker.com/r/evolinc/rmta/)
[![Docker Build Status](https://img.shields.io/docker/build/evolinc/rmta.svg)](https://hub.docker.com/r/evolinc/rmta/)
[![Docker Automated build](https://img.shields.io/docker/automated/evolinc/rmta.svg)](https://hub.docker.com/r/evolinc/rmta/)
[![Release](https://shields.beevelop.com/github/release/Evolinc/RMTA.svg?style=flat-square)](https://github.com/Evolinc/RMTA/releases)

# RMTA: Read Mapping and Transcript Assembly workflow

## Introduction

+ RMTA is a workflow that can rapidly process raw RNA-seq illumina data by mapping reads to a genome (HiSat2) and then assemble transcripts using Stringtie.
+ RMTA can process FASTq files containing paired-end or single-end reads. Alternatively, RMTA can directly process one or more sequence read archives (SRA) from NCBI using a SRA ID.
+ RMTA also supports for read alignment directly to a transcriptome using the quasi-aligner  and transcript abundance quantifier Salmon (Rob et al., 2017; Srivastava et al., 2019). Salmon maps reads to the provided transcript assembly and then counts the number of reads associated with each transcript, generating an output file (quant.sf) that can immediately be used for differential expression. **Note:** The utilization of Salmon is only appropriate when the user is wanting to rapidly test for differential expression and cannot facilitate the identification of novel genes or data visualization in a genome browser. 

Genome-guided RMTA minimally requires the following input data:

1. Reference Genome (FASTA) or Hisat2 Indexed Reference Genome (in a subdirectory)
2. RNA-Seq reads (FASTQ) - Single end or Paired end or NCBI SRA id or multiple NCBI SRA id's (list in a single column text file).

Transcriptome-guided RMTA minimally requires the following input data:

1. Reference Transcriptome (GFF3/GTF/GFF)
2. RNA-Seq reads (FASTQ) - Single end or Paired end or NCBI SRA id or multiple NCBI SRA id's (list in a single column text file).

# Availability 
### Using Docker image

Since there are several dependencies (these can be seen in Dockerfile) for running RMTA on your linux or MAC OS, we highly recommend using the available Docker image for RMTA or the [Dockerfile](https://hub.docker.com/r/evolinc/rmta/~/dockerfile/) to build an image and then use the built image. Docker can be installed on any of three platform using the instructions from [Docker website](https://docs.docker.com/engine/installation/). You can also try [Play-With-Docker](http://labs.play-with-docker.com/) for running RMTA using the below instructions 

### When using Docker
Threads and memory usage are restricted by Docker (i.e., the VM handling Docker) and thus you may have to adjust memory and CPU allowances for Docker in order to see full utilization of those resources by RMTA.

```
# Pull the image from Dockerhub
docker pull evolinc/rmta:2.6.3
```

```
# See the command line help for the image
docker run evolinc/rmta:2.6.3 -h
```

| Command Line Argument | Description |
| --------------------- | ----------- |
| -g | reference genome fasta file |
| -i | reference genome index folder |
| -A | reference genome annotation |
| -l | Library type #note that this is a lower case L |
| -1 | reads_1 |
| -2 | reads_2 |
| -U | single_reads |
| -O | /path/to/bam output folder |
| -s | SRA ID # One SRA ID or multiple SRA ID's in a file | 
| -p | Number of threads |
| -5 | 5' trim #Integers only |
| -3 | 3' trim #Integers only |
| -Q | phred64 #phred33 is default |
| -m | Minimum intron length |
| -M | Maximum intron length |
| -f | Coverage #Read per base coverage filter, integers from 0-5 only |
| -k | threshold # FPKM threshold to filter, integers from 0-5 only |
| -e | Fastqc |
| -d | Remove duplicates #intended for use with the Bowtie2 option |
| -t | Hisat2 mapping |
| -b | Bowtie2 mapping #deactivates Stringtie |
| -y | Type of reads (Single end or Paired end) #denoted as "SE" or "PE", include quotes on command line |
| -u | Feature type (Default is exon) #Feature counts option; can include any feature annotation in the provided GTF/GFF |
| -a | Gene attribute (Default is gene_id) #Feature counts option; should be the starting label in column 9 of your GTF/GFF |
| -n | Strandedness (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded)) |

```
# Run rmta on the test data. The sample data can be found in the sample_data folder in this repository (https://github.com/Evolinc/RMTA/tree/master/sample_data_arabi) 
git clone https://github.com/Evolinc/RMTA.git
cd RMTA/sample_data_arabi
```

```
# RMTA with Stringtie assembler with two paired-end fastq files with FASTqc enabled
docker run --rm -v $PWD:/data -w /data evolinc/rmta:2.6.3 -g \
genome_chr1.fa -A annotation_chr1.gtf -l "US" -n 0 -y "PE" -1 \
SRR2037320_R1.fastq.gz -1 SRR2932454_R1.fastq.gz -2 \
SRR2037320_R2.fastq.gz -2 SRR2932454_R2.fastq.gz -O final_out -p 6 \
-5 0 -3 0 -m 20 -M 50000 -t -e -u "exon" -a "gene_id" -n 0 -f 1 -k 1
```

```
# RMTA with Stringtie assembler with two single-end fastq files with no FASTqc
docker run --rm -v $PWD:/data -w /data evolinc/rmta:2.6.3 -g \
genome_chr1.fa -A annotation_chr1.gtf -l "US" -y "SE" -U SRR3464102.fastq.gz -U \
SRR3464103.fastq.gz -O final_out -p 6 -5 0 -3 0 -m 20 -M 50000 -Q -t \
-u "exon" -a "gene_id" -n 0 -f 1 -k 1
```

```
# RMTA with Stringtie assembler with one single-end fastq file
docker run --rm -v $PWD:/data -w /data evolinc/rmta:2.6.3 -g \
genome_chr1.fa -A annotation_chr1.gtf -l "US" -y "SE" -U SRR3464102.fastq.gz \
-O final_out -p 3 -5 0 -3 0 -m 20 -M 50000 -t -e -u "exon" -a "gene_id" -n 0 -f 1 -k 1
```

```
# RMTA with Stringtie assembler with one single-end fastq file with Bowtie mapper, FASTqc enabled, and duplicate reads removed
docker run --rm -v $PWD:/data -w /data evolinc/rmta:2.6.3 -g \
genome_chr1.fa -A annotation_chr1.gtf -l "US" -y "SE" -U SRR3464102.fastq.gz \
-O final_out -p 6 -5 0 -3 0 -m 20 -M 50000 -b -d -e -u "exon" -a "gene_id" -n 0 -f 1 -k 1
```

```
# One SRA id with paired-ended RNA-sequencing data running HiSat2 and FASTqc
docker run --rm -v $PWD:/data -w /data evolinc/rmta:2.6.3 -g genome_chr1.fa -A  \
annotation_chr1.gtf -l "US" -s "SRR2037320" -y "PE" \
-O final_out -p 6 -5 0 -3 0 -m 20 -M 50000 -t -e -u "exon" -a "gene_id" -n 0 -f 1 -k 1
```

```
# Multiple PE SRA's with strand-specific (forward) reads without FASTqc option
docker run --rm -v $PWD:/data -w /data evolinc/rmta:2.6.3 -g genome_chr1.fa -A annotation_chr1.gtf -l "FR" -y "PE" -s sra_id_pe.txt \
-O final_out -p 6 -5 0 -3 0 -m 20 -M 50000 -t -f 1 -k 1 -u "exon" -a "gene_id" -n 1
```

```
# transcriptome guided assembly using Salmon with FASTqc option
docker run --rm -v $PWD:/data -w /data evolinc/rmta:2.6.3 -r athal.fa.gz -s SRR2037320 -y "PE" -e -O RMTA_out -p 6
```

### Using RMTA on a High-Performance Computer (HPC)

Singularity is an open source container platform designed to be simple, fast, and secure. Singularity is optimized for EPC and HPC workloads, allowing untrusted users to run untrusted containers in a trusted way. Singularity is also good friends with Docker, hence the RMTA docker containers can be deployed on HPC using Singularity. See more documentation [here](https://sylabs.io/guides/3.4/user-guide/).

### When using Singularity

When using Singularity, you just have to replace `docker run --rm $PWD:/data -w /data evolinc/rmta:2.6.3` with `singularity run docker://evolinc/rmta:2.6.3`. For example:

```
# transcriptome guided assembly using Salmon
singularity run docker://evolinc/rmta:2.6.3 -r athal.fa.gz -s SRR2037320 -y "PE" -e -O final_out -p 6
```

### Using CyVerse Discovery Environment

The RMTA v2.6.2 app (Search for "RMTA" and then select the 2.6.2 version) is currently integrated in CyVerse’s Discovery Environment (DE) and is free to use by researchers. The complete tutorial is available at this [CyVerse wiki](https://wiki.cyverse.org/wiki/display/DEapps/RMTA+v2.6.2). CyVerse's DE is a free and easy to use GUI that simplifies many aspects of running bioinformatics analyses. If you do not currently have access to a high performance computing cluster, consider taking advantange of the DE.

# Issues
If you experience any issues with running RMTA (DE app or source code or Docker image), please open an issue on this github repo. Alternatively you can post your queries and feature requests in this [google groups](https://groups.google.com/forum/#!forum/evolinc)

# Copyright free
The sources in this [Github](https://github.com/Evolinc/RMTA) repository, are copyright free. Thus you are allowed to use these sources in which ever way you like. Here is the full [MIT](https://choosealicense.com/licenses/mit/#) license.
