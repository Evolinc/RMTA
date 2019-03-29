#!/bin/bash

# set -x
# set -e

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome>  -i <Index_folder> -A <reference_annotation> -l lib_type {-1 <left_reads> -2 <right_reads> | -U <single_reads> | -s <sra_id>} -O <output_folder for Bam files> -p num_threads -5 <integer> -3 <integer> {-q phred_33 -Q phred_64} -m min_intron -M max_intron -f <integer>"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########

  -g <reference genome fasta file>

  -i <reference genomeindex folder>

  -A <reference genome annotation>

  -l Library type

  -1 <reads_1>
               # Make sure both the paired end reads are present
  -2 <reads_2>

  -U <single_reads> # Don not use Single Reads along with Paired end reads

  -O </path/to/bam output folder>

  -s SRA ID # One SRA ID or multiple SRA ID's in a file 

  -p Number of threads
  
  -5 5' trim 
 
  -3 3' trim

  -q phred33

  -Q phred64 

  -m Minimum intron length

  -M Maximum intron length

  -f threshold # FPKM threshold to filter

EOF
    exit 0
}

quality_64=0
fastqc=0
hisat=0
bowtie=0
duplicates=0

while getopts ":hg:i:A:l:1:2:U:O:s:p:5:3:f:Qedtbm:M:k:" opt; do
  case $opt in
    g)
    referencegenome=$OPTARG # Reference genome file
     ;;
    i)
    index_folder=$OPTARG # Input folder
     ;;
    A)
    referenceannotation=$OPTARG # Reference genome annotation
     ;;
    l)
     lib_type=$OPTARG # Library type
     ;;
    1)
    left_reads+=("$OPTARG") # Left reads
     ;;
    2)
    right_reads=("$OPTARG") # Right reads
     ;;
    U)
    single_reads+=("$OPTARG") # single end reads
     ;;
    O)
    bam_out=$OPTARG # Samoutput file
     ;;
    s)
    sra_id=$OPTARG # SRA ID or SRA ID's in a file
     ;;
    p)
    num_threads=$OPTARG # Number of threads
     ;;
    5)
    five_trim=$OPTARG # 5' trim
     ;;
    3)
    three_trim=$OPTARG # 3' trim
     ;;
    Q)
    quality_64=$OPTARG # Phread 64
     ;;
    m)
    min_intl=$OPTARG # Minimum intron length
     ;;
    M)
    max_intl=$OPTARG # Maximum intron length
     ;;
    f)
    threshold=$OPTARG # Coverage/base filter that you would like to apply to identified transcripts
     ;;
    k)
    fthreshold=$OPTARG # FPKM filter that you would like to apply to identified transcripts
     ;;
    e)
    fastqc=$OPTARG # Fastqc
     ;;
    d)
    duplicates=$OPTARG # Duplicates removal
     ;;
    t)
    hisat=$OPTARG # Hisat2
     ;;
    b)
    bowtie=$OPTARG # Bowtie2
     ;;
    h)
    usage
     exit 1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# ############################################################################################################################################################################################################################
## Functions ###
# ############################################################################################################################################################################################################################

#Define parameters for coverage/base filtering
param5=5
param4=4
param3=3
param2=2
param1=1
param0=0

##################
### No SRA #######
##################

### Paired end #####

# Coverage cut-off function for paired end reads

coverge_cuffoff_non_SRA()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "4\.' -e 'cov "3\.' -e 'cov "2\.' -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "3\.' -e 'cov "2\.' -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "2\.' -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "1\.' -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'cov "0\.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# FPKM cut-off function for paired end reads

fpkm_cuffoff_non_SRA()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "4\.' -e 'FPKM "3\.' -e 'FPKM "2\.' -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "3\.' -e 'FPKM "2\.' -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "2\.' -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "1\.' -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "0\.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename3.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename3.gtf | grep -e 'FPKM "0\.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf >"$bam_out"/$filename3.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Stringtie function for paired end reads

stringtie_non_SRA() 
{
    if [ "$hisat" != 0 ]; then
      echo "###################"
      echo "Running Stringtie"
      echo "####################"
      echo "stringtie -G $referenceannotation $filename3.sorted.bam -o $filename3.gtf -p $num_threads"
      stringtie -G $referenceannotation $filename3.sorted.bam -o $filename3.gtf -p $num_threads
      mv $filename3.sorted.bam $filename3.sorted.bam.bai $filename3.gtf "$bam_out"
      coverge_cuffoff_non_SRA
      fpkm_cuffoff_non_SRA 
      var=$(grep -vFf listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf | wc -l)
      if [ "$var" -eq 0 ]; then
        rm listtoremove.txt "$bam_out"/$filename3.gtf.filtered.gtf
      else
        echo "###################"
        echo "Running Cuffcompare"
        echo "####################"
        echo "cuffcompare "$bam_out"/$filename3.gtf.filtered.gtf -r $referenceannotation -o $filename3"
        cuffcompare "$bam_out"/$filename3.gtf.filtered.gtf -r $referenceannotation -o $filename3
        rm listtoremove.txt
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
      fi
    fi
}

### single end ####

# Coverage cut-off function for single end reads

coverge_cuffoff_non_SRA_single()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}
 
# FPKM cut-off function for paired end reads

fpkm_cuffoff_non_SRA_single()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "4.' -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $filename.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$filename.gtf | grep -e 'FPKM "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$filename.gtf >"$bam_out"/$filename.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Stringtie function for single end reads

stringtie_non_SRA_single() 
{
    if [ "$hisat" != 0 ]; then
      echo "###################"
      echo "Running Stringtie"
      echo "####################"
      echo "stringtie -G $referenceannotation $filename.sorted.bam -o $filename.gtf -p $num_threads"
      stringtie -G $referenceannotation $filename.sorted.bam -o $filename.gtf -p $num_threads
      mv $filename.sorted.bam $filename.sorted.bam.bai $filename.gtf "$bam_out"
      coverge_cuffoff_non_SRA_single
      fpkm_cuffoff_non_SRA_single
      var=$(grep -vFf listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf | wc -l)
      if [ "$var" -eq 0 ]; then
        rm listtoremove.txt "$bam_out"/$filename.gtf.filtered.gtf
      else
        echo "###################"
        echo "Running Cuffcompare"
        echo "####################"
        echo "cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename"
        cuffcompare "$bam_out"/$filename.gtf.filtered.gtf -r $referenceannotation -o $filename
        rm listtoremove.txt
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
      fi
    fi
}

##################
### Single SRA ###
##################

# Coverage cut-off function for single SRAs

coverge_cuffoff_SRA_single()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# FPKM cut-off cunction for Single SRAs

fpkm_cuffoff_SRA_single()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "4.' -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $sra_id.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi             
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$sra_id.gtf | grep -e 'FPKM "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf >"$bam_out"/$sra_id.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Stringtie function for single SRAs

stringtie_SRA_single()
{
    if [ "$hisat" != 0 ]; then
      echo "###################"
      echo "Running Stringtie"
      echo "####################"
      echo "stringtie -G $referenceannotation $sra_id.sorted.bam -o $sra_id.gtf -p $num_threads"
      stringtie -G $referenceannotation $sra_id.sorted.bam -o $sra_id.gtf -p $num_threads
      mv $sra_id.sorted.bam $sra_id.sorted.bam.bai $sra_id.gtf "$bam_out"
      coverge_cuffoff_SRA_single
      fpkm_cuffoff_SRA_single  
      var=$(grep -vFf listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf | wc -l)
      if [ "$var" -eq 0 ]; then
        rm listtoremove.txt "$bam_out"/$sra_id.gtf.filtered.gtf
      else
        echo "###################"
        echo "Running Cuffcompare"
        echo "####################"
        echo "cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id"
        cuffcompare "$bam_out"/$sra_id.gtf.filtered.gtf -r $referenceannotation -o $sra_id
        rm listtoremove.txt
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
      fi
    fi
}

##################
### Multi SRA ###
##################

# Coverage cut-off function for multiple SRAs

coverge_cuffoff_SRA_multi()
{
    if [ "$threshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "4.' -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$threshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "3.' -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "2.' -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$threshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "1.' -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$threshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# FPKM cut-off function for multiple SRAs

fpkm_cuffoff_SRA_multi()
{
    if [ "$fthreshold" -eq "$param5" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "4.' -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param5" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$fthreshold" -eq "$param4" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "3.' -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param4" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param3" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "2.' -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param3" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi      
    elif [ "$fthreshold" -eq "$param2" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "1.' -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param2" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param1" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'FPKM "0.' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
       var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
       if [ "$var" -eq 0 ]; then
          echo No transcripts have FPKM values exceeding your "$param1" cut-off have been found in $f.gtf.filtered.gtf. Try lowering your cut-off 1>&2
       fi       
    elif [ "$fthreshold" -eq "$param0" ]; then
       grep " transcript" "$bam_out"/$f.gtf | grep -e 'cov "0.000' | cut -f 9 | cut -d " " -f 4 | sort -u > listtoremove.txt
       grep -vFf listtoremove.txt "$bam_out"/$f.gtf >"$bam_out"/$f.gtf.filtered.gtf
    else
        echo "Invalid coverage parameter. Please select a whole number between 0-5"
        exit
    fi
}

# Success message

success_message()
{
  if [ ! -z "$index_folder" ]; then
    rm $fbname*
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    mkdir index
    mv $fbname* index
    mv mapped.txt metrics.txt "$bam_out"
  fi
  if [ "$hisat" != 0 ]; then
    mkdir Feature_counts && mv feature_counts.* Feature_counts
    mv Feature_counts "$bam_out"
  fi
  echo "##############################"
  echo "Pipeline executed successfully"
  echo "##############################"
}

# Feature counts

featurecounts()
{
    if [ "$lib_type" == "F" ] || [ "$lib_type" == "R" ]; then
      featureCounts -T $num_threads -t exon -g gene_id -a $referenceannotation -o feature_counts.txt "$bam_out"/*.sorted.bam
    elif [ "$lib_type" == "FR" ] || [ "$lib_type" == "RF" ]; then
      featureCounts -p -T $num_threads -t exon -g gene_id -a $referenceannotation -o feature_counts.txt "$bam_out"/*.sorted.bam
    fi
    awk '{$2=$3=$4=$5=$6=""; print $0}' OFS='\t' feature_counts.txt | grep -v "#" | sed 's/\t\+/\t/g;s/^\t//' > temp.txt && mv temp.txt feature_counts.txt
}

# Duplicate removal for paired end reads

duplicates_paired()
{
    echo "sambamba sort --tmpdir=temp -t $num_threads -o $filename3.sorted.bam $filename3.bam"
    sambamba sort --tmpdir=temp -t $num_threads -o $filename3.sorted.bam $filename3.bam
    rm -r temp $filename3.bam
    if [ "$duplicates" != 0 ]; then
      rmdup=$(basename $filename3.sorted.bam ".sorted.bam")
      picard MarkDuplicates I=$filename3.sorted.bam O=$rmdup."sorted.rmdup.bam" ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
      mv $filename3.sorted.bam $filename3.sorted.bam.bai $rmdup."sorted.rmdup.bam" "$bam_out"
    fi
}

# Stringtie function for multiple SRAs

stringtie_SRA_multi()
{
    if [ "$hisat" != 0 ]; then
      echo "###################"
      echo "Running Stringtie"
      echo "####################"
      echo "stringtie -G $referenceannotation $f.sorted.bam -o $f.gtf -p $num_threads"
      stringtie -G $referenceannotation $f.sorted.bam -o $f.gtf -p $num_threads
      mv $f.sorted.bam $f.sorted.bam.bai $f.gtf "$bam_out"
      coverge_cuffoff_SRA_multi
      fpkm_cuffoff_SRA_multi
      var=$(grep -vFf listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf | wc -l)
      if [ "$var" -eq 0 ]; then
        rm listtoremove.txt "$bam_out"/$f.gtf.filtered.gtf
      else
        echo "###################"
        echo "Running Cuffcompare"
        echo "####################"
        echo "cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f"
        cuffcompare "$bam_out"/$f.gtf.filtered.gtf -r $referenceannotation -o $f
        rm listtoremove.txt
        mv *.tracking *.loci *.combined.gtf *.stats "$bam_out"
      fi
    fi
}

# Fq.gz and Fastq_gz for paired end reads

paired_fq_gz()
{
    filename=$(basename "$f" ".fq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    fi    
    echo "samtools view -Sbo $filename3.bam $filename3.sam"
    samtools view -Sbo $filename3.bam $filename3.sam
    duplicates_paired
    stringtie_non_SRA
    if [ "$hisat" != 0 ]; then
      featurecounts
    fi    
    rm $filename3.sam
}

paired_fq()
{
    filename=$(basename "$f" ".fq")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    fi
    echo "samtools view -Sbo $filename3.bam $filename3.sam"
    samtools view -Sbo $filename3.bam $filename3.sam
    duplicates_paired
    stringtie_non_SRA
    if [ "$hisat" != 0 ]; then
      featurecounts
    fi
    rm $filename3.sam
}

paired_fastq_gz()
{
    filename=$(basename "$f" ".fastq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    fi    
    echo "samtools view -Sbo $filename3.bam $filename3.sam"
    samtools view -Sbo $filename3.bam $filename3.sam
    duplicates_paired
    stringtie_non_SRA
    if [ "$hisat" != 0 ]; then
      featurecounts
    fi
    rm $filename3.sam
}

paired_fastq()
{
    filename=$(basename "$f" ".fastq")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
    fi
    echo "samtools view -Sbo $filename3.bam $filename3.sam"
    samtools view -Sbo $filename3.bam $filename3.sam
    duplicates_paired
    stringtie_non_SRA
    if [ "$hisat" != 0 ]; then
      featurecounts
    fi    
    rm $filename3.sam
}

single_end()
{
    extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
    filename=$(basename "$f" ".$extension")
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $filename "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
        echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
    fi
    echo "samtools view -Sbo $filename.bam $filename.sam"
    samtools view -Sbo $filename.bam $filename.sam
    echo "sambamba sort --tmpdir=temp -t $num_threads -o $filename.sorted.bam $filename.bam"
    sambamba sort --tmpdir=temp -t $num_threads -o $filename.sorted.bam $filename.bam
    rm -r temp $filename.bam
    if [ "$duplicates" != 0 ]; then
      rmdup=$(basename $filename.sorted.bam ".sorted.bam")
      picard MarkDuplicates I=$filename.sorted.bam O=$rmdup."sorted.rmdup.bam" ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
      mv $filename.sorted.bam $filename.sorted.bam.bai $rmdup".sorted.rmdup.bam" "$bam_out"
    fi
    stringtie_non_SRA_single
    if [ "$hisat" != 0 ]; then
      featurecounts
    fi
    rm $filename.sam
}

# ############################################################################################################################################################################################################################
# # Reference genome/Index
# ############################################################################################################################################################################################################################

if [ "$hisat" != 0 ] && [ "$bowtie" == 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
        cp $i .
        fbname=$(basename "$i" .ht2 | cut -d. -f1)
    done
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    echo "####################################"
    echo "Starting Reference genome build"
    echo "####################################"
    echo "hisat2-build $referencegenome ref_genome -p $num_threads"
    hisat2-build $referencegenome ref_genome -p $num_threads
    echo "fbname=$(basename "ref_genome" .ht2 | cut -d. -f1)"
    fbname=$(basename "ref_genome" .ht2 | cut -d. -f1)
  fi
elif [ "$hisat" == 0 ] && [ "$bowtie" != 0 ]; then
  if [ ! -z "$index_folder" ]; then
    for i in $index_folder/*; do
        cp $i .
        fbname=$(basename "$i" .ht2 | cut -d. -f1)
    done
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    echo "####################################"
    echo "Starting Reference genome build"
    echo "####################################"
    echo "bowtie2-build $referencegenome ref_genome --threads $num_threads"
    bowtie2-build $referencegenome ref_genome -p $num_threads
    echo "fbname=$(basename "ref_genome" .ht2 | cut -d. -f1)"
    fbname=$(basename "ref_genome" .bt2 | cut -d. -f1)
  fi
fi

# ############################################################################################################################################################################################################################
# # Mapping
# ############################################################################################################################################################################################################################

# Paired end reads

if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ]; then
    mkdir "$bam_out"
    numb=$(ls "${left_reads[@]}" | wc -l)
    for f in "${left_reads[@]}"; do
      extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
      if [[ "$extension" =~ "fq.gz" ]]; then
        paired_fq_gz
      elif [[ "$extension" =~ "fastq.gz" ]]; then
        paired_fastq_gz
      elif [[ "$extension" =~ "fq" ]]; then
        echo "gzip" "$f"
        paired_fq
      elif [[ "$extension" =~ "fastq" ]]; then
        echo "gzip" "$f"
        paired_fastq
      elif [ "$extension" != "fastq" ] || [ "$extension" != "fq" ] || [ "$extension" != "fastq.gz" ] || [ "$extension" != "fq.gz" ]; then
        echo "The extension" "$extension" "is not supported. Only .fq, .fq.gz, .fastq, .fastq.gz are only supported" 1>&2        
        exit 64
      fi 
    done
    success_message

# Single end reads

elif [ ! -z "$single_reads" ]; then
    mkdir "$bam_out"
    numb=$(ls "${single_reads[@]}" | wc -l)
    for f in "${single_reads[@]}"; do
      single_end
    done
    success_message

      
# SRA

elif [ ! -z $sra_id ]; then
  if [[ -f $sra_id ]]; then
    if [ "$fastqc" != 0 ]; then
      mkdir -p $bam_out/Fasqc_out
    else
      mkdir "$bam_out" 
    fi  
    while read f; do
      echo "------------------------------------------" >> mapped.txt
      echo "### Mapping percentages of" $f "###" >> mapped.txt
      echo "------------------------------------------" >> mapped.txt
      echo "#######################"
      echo "Running Hisat2 mapping"
      echo "#######################"
      if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then #(-Q -h)
          echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
          hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then #(-Q -b)
          echo "bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
          bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default #(-h)
          echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
          hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then #(-b)
          echo "bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
          bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt
      fi
      samtools view -Sbo $f.bam $f.sam
      rm $f.sam
      echo "sambamba sort --tmpdir=temp -t $num_threads -o $f.sorted.bam $f.bam"
      sambamba sort --tmpdir=temp -t $num_threads -o $f.sorted.bam $f.bam
      if [ "$duplicates" != 0 ]; then
        rmdup=$(basename $f.sorted.bam ".sorted.bam")
        picard MarkDuplicates I=$f.sorted.bam O=$rmdup."sorted.rmdup.bam" ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
        mv $f.sorted.bam $f.sorted.bam.bai $rmdup."sorted.rmdup.bam" "$bam_out"
      fi
      stringtie_SRA_multi
      if [ "$fastqc" != 0 ]; then
        if [ "$lib_type" == "F" ] || [ "$lib_type" == "R" ]; then
          picard SamToFastq I="$f".bam FASTQ="$f".fastq
          fastqc "$f".fastq
          mkdir "$f".fastqc_out
          mv "$f".fastq *.zip *.html "$f".fastqc_out
        elif [ "$lib_type" == "FR" ] || [ "$lib_type" == "RF" ]; then
          picard SamToFastq I="$f".bam FASTQ="$f"."R1".fastq SECOND_END_FASTQ="$f"."R2".fastq
          fastqc "$f"."R1".fastq "$f"."R2".fastq
          mkdir "$f".fastqc_out
          mv "$f".*.fastq *.zip *.html "$f".fastqc_out
        fi
      fi
    done < "$sra_id"
      if [ "$fastqc" != 0 ]; then
        mv *.fastqc_out "$bam_out/Fasqc_out"
      fi
      if [ "$hisat" != 0 ]; then
        featurecounts
      fi
      rm -r temp *.bam
      success_message
  else    
    mkdir "$bam_out"
    echo "------------------------------------------" >> mapped.txt
    echo "### Mapping percentages of" $sra_id "###" >> mapped.txt
    echo "------------------------------------------" >> mapped.txt
    echo "#######################"
    echo "Running Hisat2 mapping"
    echo "#######################"
    if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
        echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
        echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
        hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
        echo "bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
        bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
    fi
    samtools view -Sbo $sra_id.bam $sra_id.sam
    rm $sra_id.sam
    echo "sambamba sort --tmpdir=temp -t $num_threads -o $sra_id.sorted.bam $sra_id.bam"
    sambamba sort --tmpdir=temp -t $num_threads -o $sra_id.sorted.bam $sra_id.bam
    if [ "$duplicates" != 0 ]; then
      rmdup=$(basename $sra_id.sorted.bam ".sorted.bam")
      picard MarkDuplicates I=$sra_id.sorted.bam O=$rmdup."sorted.rmdup.bam" ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
      mv $sra_id.sorted.bam $sra_id.sorted.bam.bai $rmdup."sorted.rmdup.bam" "$bam_out"
    fi
    stringtie_SRA_single
    if [ "$fastqc" != 0 ]; then
      if [ "$lib_type" == "F" ] || [ "$lib_type" == "R" ]; then
        picard SamToFastq I=$sra_id.bam FASTQ="$sra_id".fastq
        fastqc "$sra_id".fastq
        mkdir fastqc_out
        mv "$sra_id".fastq *.zip *.html fastqc_out
        mv fastqc_out "$bam_out/Fasqc_out"
      elif [ "$lib_type" == "FR" ] || [ "$lib_type" == "RF" ]; then
        picard SamToFastq I=$sra_id.bam FASTQ="$sra_id"."R1".fastq SECOND_END_FASTQ="$sra_id"."R2".fastq
        echo "fastqc "$sra_id"."R1".fastq "$sra_id"."R2".fastq"
        fastqc "$sra_id"."R1".fastq "$sra_id"."R2".fastq
        mkdir fastqc_out
        mv "$sra_id".*.fastq *.zip *.html fastqc_out
        mv fastqc_out "$bam_out/Fasqc_out"
      fi
    fi
    if [ "$hisat" != 0 ]; then
      featurecounts
    fi
    rm -r temp $sra_id.bam
    success_message
  fi
fi
##### End ########