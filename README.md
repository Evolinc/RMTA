# RMTA-1.0

RMTA (Read Mapping and Transcript Assembly) is a workflow that can process raw RNA-seq data either paired-end or single-end or from NCBI’s sequence read archive (SRA) directly and generate transcript assembly

RMTA minimally requires the following input data

1. Reference Genome (FASTA) or Hisat2 Indexed Reference Genome (in a folder)
2. Reference Transcriptome (GFF3)
3. RNA-Seq reads (FASTQ) - Single end or Paired end or NCBI SRA id or multiple NCBI SRA id's

# Availability 
### Using Docker image

Since there are several dependencies (these can be seen in Dockerfile) for running Evolinc-I on your linux or MAC OS, we highly recommend using the available Docker image for RMTA or the [Dockerfile](https://hub.docker.com/r/evolinc/rmta/~/dockerfile/) to build an image and then use the built image.

```
# Pull the image from Dockerhub
docker pull evolinc/rmta:1.0

# See the command line help for the image
docker run evolinc/rmta:1.0 -h

# Run rmta on the test data. The sample data can be found in the sample_data folder in this (https://github.com/Evolinc/RMTA/tree/master/sample_data) 
git clone https://github.com/Evolinc/RMTA.git
cd RMTA

# Hisat2 + Stringtie with two fastq files
docker run --rm -v $PWD:/data -w /data evolinc/rmta:1.0 -g \
Sorghum_bicolor.Sorbi1.20.dna.toplevel_chr8.fa -A \
Sorghum_bicolor.Sorbi1.20_chr8.gtf -l "FR" -1 sample_1_R1.fq.gz -1 \
sample_2_R1.fq.gz -2 sample_1_R2.fq.gz -2 sample_2_R2.fq.gz -O final_out \
-p 6 -5 0 -3 0 -m 20 -M 50000 -q -t -f 2 

# Hisat2 + Cufflinks with two fastq files
docker run --rm -v $PWD:/data -w /data evolinc/rmta:1.0 -g \
Sorghum_bicolor.Sorbi1.20.dna.toplevel_chr8.fa -A \
Sorghum_bicolor.Sorbi1.20_chr8.gtf -l "FR" -1 sample_1_R1.fq.gz -1 \
sample_2_R1.fq.gz -2 sample_1_R2.fq.gz -2 sample_2_R2.fq.gz -O final_out \
-p 6 -5 0 -3 0 -m 20 -M 50000 -q -t -f 2 

# One SRA id
docker run --rm -v $PWD:/data -w /data evolinc/rmta:1.0 -g Sorghum_bicolor.Sorbi1.20.dna.toplevel_chr8.fa -A Sorghum_bicolor.Sorbi1.20_chr8.gtf -l "FR" -s "SRR3993757" \
-O final_out -p 6 -5 0 -3 0 -m 20 -M 50000 -q -t -f 2 

# Multiple SRA's
docker run --rm -v $PWD:/data -w /data evolinc/rmta:1.0 -g Sorghum_bicolor.Sorbi1.20.dna.toplevel_chr8.fa -A Sorghum_bicolor.Sorbi1.20_chr8.gtf -l "FR" -s sra_id.txt \
-O final_out -p 6 -5 0 -3 0 -m 20 -M 50000 -q -t -f 2 
```

# Issues
If you experience any issues with running Evolinc-I (DE app or source code or Docker image), please open an issue on this github repo.

# Copyright free
The sources in this [Github](https://github.com/Evolinc/Evolinc-I) repository, are copyright free. Thus you are allowed to use these sources in which ever way you like. Here is the full [MIT](https://choosealicense.com/licenses/mit/#) license.
