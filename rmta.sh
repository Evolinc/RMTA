#!/bin/bash

# set -x
# set -e

usage() {
      echo ""
      echo "Usage : sh $0 -g <reference_genome>  -i <Index_folder> -A <reference_annotation> -l lib_type {-1 <left_reads> -2 <right_reads> | -U <single_reads> | -s <sra_id>} -O <output_folder for Bam files> -p num_threads -5 <integer> -3 <integer> -Q phred_64 (default 33) -m min_intron -M max_intron -f <integer> -k <integer> {-e fastqc} -h hisat2 -b bowtie2 {-d duplicates}"
      echo ""

cat <<'EOF'
  
  ###### Command line options ##########

  -g <reference genome fasta file>

  -i <reference genomeindex folder>

  -A <reference genome annotation>

  -l Library type #note that this is a lower case L

  -1 <reads_1>
               # Ends with R1 and is in the same order as reverse reads
  -2 <reads_2>
               # Ends with R2, must be present, and is in the same order as forward reads
  -U <single_reads> # Do not use Single Reads along with Paired end reads

  -O </path/to/bam output folder>

  -s SRA ID # One SRA ID or multiple SRA ID's in a file 

  -p Number of threads
  
  -5 5' trim #Integers only
 
  -3 3' trim #Integers only

  -Q phred64 #phred33 is default if this flag is not activated

  -m Minimum intron length

  -M Maximum intron length

  -f Coverage #Read per base coverage filter, integers from 0-5 only

  -k threshold # FPKM threshold to filter, integers from 0-5 only

  -e Fastqc

  -d Remove duplicates #intended for use with the Bowtie2 option

  -t Hisat2 mapping

  -b Bowtie2 mapping #deactivates Stringtie

  -y Type of reads (Single end or Paired end) #denoted as "SE" or "PE", include quotes on command line

  -u Feature type (Default is exon) #Feature counts option; can include any feature annotation in the provided GTF/GFF

  -a Gene attribute (Default is gene_id) #Feature counts option; should be the starting label in column 9 of your GTF/GFF

  -n Strandedness (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded)) #Feature counts option; separate strand information for Featurecounts

EOF
    exit 0
}

quality_64=0
fastqc=0
hisat=0
bowtie=0
duplicates=0
referenceannotation=0
transcriptome=0
referencegenome=0

while getopts ":hg:i:A:r:l:1:2:U:O:s:p:5:3:f:QAedtbm:M:k:y:u:a:n:" opt; do
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
    r)
    transcriptome=$OPTARG # transcriptome
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
    y)
    seq_type=$OPTARG # Type of Sequence data (SE or PE. Mainly needed for SRA and featurecounts)
     ;;
    u) 
    feature_type=$OPTARG # Feature type (Default is exon)
     ;;
    a) 
    gene_attribute=$OPTARG # (Default is gene_id)
     ;;
    n)
    strandedness=$OPTARG # (Default is 0 (unstranded), 1 (stranded), 2 (reversely stranded))
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
### Functions ###
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
      if [ "$referenceannotation" != 0 ]; then
        echo "stringtie -G $referenceannotation "$bam_out"/$filename3.sorted.bam -o "$bam_out"/$filename3.gtf -p $num_threads"
        stringtie -G $referenceannotation "$bam_out"/$filename3.sorted.bam -o "$bam_out"/$filename3.gtf -p $num_threads
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
          rm listtoremove.txt *.tracking *.loci *.stats 
          mv *.combined.gtf "$bam_out"
          rm "$bam_out"/*refmap "$bam_out"/*tmap
        fi
      else
        echo "stringtie "$bam_out"/$filename3.sorted.bam -o "$bam_out"/$filename3.gtf -p $num_threads"
        stringtie "$bam_out"/$filename3.sorted.bam -o "$bam_out"/$filename3.gtf -p $num_threads
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
      if [ "$referenceannotation" != 0 ]; then
        echo "stringtie -G $referenceannotation $filename.sorted.bam -o $filename.gtf -p $num_threads"
        stringtie -G $referenceannotation "$bam_out"/$filename.sorted.bam -o "$bam_out"/$filename.gtf -p $num_threads
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
          rm listtoremove.txt *.tracking *.loci *.stats 
          mv *.combined.gtf "$bam_out"
          rm "$bam_out"/*refmap "$bam_out"/*tmap
        fi
      else
        echo "stringtie $filename.sorted.bam -o $filename.gtf -p $num_threads"
        stringtie "$bam_out"/$filename.sorted.bam -o "$bam_out"/$filename.gtf -p $num_threads
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

# FPKM cut-off function for Single SRAs

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
      if [ "$referenceannotation" != 0 ]; then
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
          rm listtoremove.txt *.tracking *.loci *.stats 
          mv *.combined.gtf "$bam_out"
          rm "$bam_out"/*refmap "$bam_out"/*tmap
        fi
      else
        echo "stringtie $sra_id.sorted.bam -o $sra_id.gtf -p $num_threads"
        stringtie $sra_id.sorted.bam -o $sra_id.gtf -p $num_threads
        mv $sra_id.sorted.bam $sra_id.sorted.bam.bai $sra_id.gtf "$bam_out"
      fi
    else 
      echo "###################"
      echo "Bowtie was selected, skipping Stringtie"
      echo "####################"
      mv $sra_id.sorted.bam $sra_id.sorted.bam.bai "$bam_out"
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

# Success message #

success_message()
{
  if [ ! -z "$index_folder" ]; then
    rm $fbname*
  elif [ ! -z "$referencegenome" ] && [ -z "$index_folder" ]; then
    mkdir index
    mv $fbname* index
    mv mapped.txt metrics.txt "$bam_out"
  fi
  if [ "$referenceannotation" != 0 ]; then
    if [ "$hisat" != 0 ]; then
      mkdir Feature_counts && mv feature_counts.* Feature_counts
      mv Feature_counts "$bam_out"
    fi
  fi
  echo "##############################"
  echo "Pipeline executed successfully"
  echo "##############################"
}

# Feature counts

featurecounts()
{
    if [ "$referenceannotation" != 0 ]; then
      if [ "$seq_type" == "SE" ]; then
        echo "featureCounts -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o feature_counts.txt "$bam_out"/*.sorted.bam"
        featureCounts -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o feature_counts.txt "$bam_out"/*.sorted.bam
      elif [ "$seq_type" == "PE" ]; then
        echo "featureCounts -p -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o feature_counts.txt "$bam_out"/*.sorted.bam"
        featureCounts -p -T $num_threads -t $feature_type -g $gene_attribute -s $strandedness -a $referenceannotation -o feature_counts.txt "$bam_out"/*.sorted.bam
      fi
      awk '{$2=$3=$4=$5=$6=""; print $0}' OFS='\t' feature_counts.txt | grep -v "#" | sed 's/\t\+/\t/g;s/^\t//' > temp.txt && mv temp.txt feature_counts.txt
      sed -i "s/$bam_out\///g" feature_counts.txt
    fi
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
    else
      mv $filename3.sorted.bam $filename3.sorted.bam.bai "$bam_out"
    fi
}

# Stringtie function for multiple SRAs

stringtie_SRA_multi()
{
    if [ "$hisat" != 0 ]; then
      echo "###################"
      echo "Running Stringtie"
      echo "####################"
      if [ "$referenceannotation" != 0 ]; then
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
          rm listtoremove.txt *.tracking *.loci *.stats 
          mv *.combined.gtf "$bam_out"
          rm "$bam_out"/*refmap "$bam_out"/*tmap
        fi
      else
        echo "stringtie $f.sorted.bam -o $f.gtf -p $num_threads"
        stringtie $f.sorted.bam -o $f.gtf -p $num_threads
        mv $f.sorted.bam $f.sorted.bam.bai $f.gtf "$bam_out"
      fi
    else 
      echo "###################"
      echo "Bowtie was selected, skipping Stringtie"
      echo "####################"
      mv $sra_id.sorted.bam $sra_id.sorted.bam.bai "$bam_out"
    fi
}

# Fq.gz and Fastq_gz for paired end reads

paired_fq_gz()
{
    filename=$(basename "$f" ".fq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    if [ "$referencegenome" != 0 ] || [ "$index_folder" != 0 ]; then
      if [ "$transcriptome" == 0 ]; then
        mkdir $bam_out
        if [ "$fastqc" != 0 ]; then
          DIRECTORY=$bam_out/Fastqc_out
          if [ ! -d "$DIRECTORY" ]; then
            mkdir $DIRECTORY
          fi
          fastqc ${filename}.fq.gz ${filename2}.fq.gz
          mkdir "$filename3".fastqc_out
          mv *.zip *.html "$filename3".fastqc_out
          mv "$filename3".fastqc_out $DIRECTORY
        fi
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Mapping step"
        echo "#######################"
        if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
        elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
        fi    
        echo "samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam"
        samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam
        duplicates_paired
        stringtie_non_SRA
        if [ "$hisat" != 0 ]; then
          featurecounts
        fi
        rm $filename3.sam
      else
        if [ "$transcriptome" != 0 ]; then
          mkdir -p $bam_out/Salmon_counts
          if [ "$fastqc" != 0 ]; then
            DIRECTORY=$bam_out/Fastqc_out
            if [ ! -d "$DIRECTORY" ]; then
              mkdir $DIRECTORY
            fi
            fastqc ${filename}.fq.gz ${filename2}.fq.gz
            mkdir "$filename3".fastqc_out
            mv *.zip *.html "$filename3".fastqc_out
            mv "$filename3".fastqc_out $DIRECTORY
          fi
          salmon quant -i salmon_index -l A -1 ${filename}.fq.gz -2 ${filename2}.fq.gz -p $num_threads --validateMappings -o $filename3.quant
          mv $filename3.quant $bam_out/Salmon_counts
        fi
      fi
    fi
}

paired_fq()
{
    filename=$(basename "$f" ".fq")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    if [ "$referencegenome" != 0 ] || [ "$index_folder" != 0 ]; then
      if [ "$transcriptome" == 0 ]; then
        if [ ! -d "$bam_out" ]; then
            mkdir $bam_out
        fi
        if [ "$fastqc" != 0 ]; then
          DIRECTORY=$bam_out/Fastqc_out
          if [ ! -d "$DIRECTORY" ]; then
            mkdir $DIRECTORY
          fi
          fastqc ${filename}.fq ${filename2}.fq
          mkdir "$filename3".fastqc_out
          mv *.zip *.html "$filename3".fastqc_out
          mv "$filename3".fastqc_out $DIRECTORY
        fi
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Mapping step"
        echo "#######################"
        if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
        elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
        fi    
        echo "samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam"
        samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam
        duplicates_paired
        stringtie_non_SRA
        if [ "$hisat" != 0 ]; then
          featurecounts
        fi
        rm $filename3.sam
      else
        if [ "$transcriptome" != 0 ]; then
          mkdir -p $bam_out/Salmon_counts
          if [ "$fastqc" != 0 ]; then
            DIRECTORY=$bam_out/Fastqc_out
            if [ ! -d "$DIRECTORY" ]; then
              mkdir $DIRECTORY
            fi
            fastqc ${filename}.fq ${filename2}.fq
            mkdir "$filename3".fastqc_out
            mv *.zip *.html "$filename3".fastqc_out
            mv "$filename3".fastqc_out $DIRECTORY
          fi
          salmon quant -i salmon_index -l A -1 ${filename}.fq -2 ${filename2}.fq -p $num_threads --validateMappings -o $filename3.quant
          mv $filename3.quant $bam_out/Salmon_counts
        fi
      fi
    fi
}

paired_fastq_gz()
{
    filename=$(basename "$f" ".fastq.gz")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    if [ "$referencegenome" != 0 ] || [ "$index_folder" != 0 ]; then
      if [ "$transcriptome" == 0 ]; then
        if [ ! -d "$bam_out" ]; then
              mkdir $bam_out
        fi
  	    if [ "$fastqc" != 0 ]; then
          DIRECTORY=$bam_out/Fastqc_out
          if [ ! -d "$DIRECTORY" ]; then
            mkdir $DIRECTORY
          fi
  	      fastqc ${filename}.fastq.gz ${filename2}.fastq.gz
  	      mkdir "$filename3".fastqc_out
  	      mv *.zip *.html "$filename3".fastqc_out
  	      mv "$filename3".fastqc_out $DIRECTORY
  	    fi
  	    echo "------------------------------------------" >> mapped.txt
  	    echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
  	    echo "------------------------------------------" >> mapped.txt
  	    echo "#######################"
  	    echo "Running Mapping step"
  	    echo "#######################"
  	    if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
  	      if [ "$lib_type" == "US" ]; then
  	        echo "hisat2 -x $fbname -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
  	        hisat2 -x $fbname -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
  	      else
  	        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
  	        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
  	      fi
  	    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
  	        echo "bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
  	        bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
  	    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
  	      if [ "$lib_type" == "US" ]; then
  	        echo "hisat2 -x $fbname -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
  	        hisat2 -x $fbname -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
  	      else
  	        echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
  	        hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
  	      fi
  	    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
  	        echo "bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
  	        bowtie2 -x ref_genome -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
  	    fi    
  	    echo "samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam"
  	    samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam
  	    duplicates_paired
  	    stringtie_non_SRA
  	    if [ "$hisat" != 0 ]; then
  	      featurecounts
  	    fi
  	    rm $filename3.sam
	    else
        if [ "$transcriptome" != 0 ]; then
          mkdir -p $bam_out/Salmon_counts
          if [ "$fastqc" != 0 ]; then
            DIRECTORY=$bam_out/Fastqc_out
            if [ ! -d "$DIRECTORY" ]; then
              mkdir $DIRECTORY
            fi
            fastqc ${filename}.fastq.gz ${filename2}.fastq.gz
            mkdir "$filename3".fastqc_out
            mv *.zip *.html "$filename3".fastqc_out
            mv "$filename3".fastqc_out $DIRECTORY
          fi
      	  salmon quant -i salmon_index -l A -1 ${filename}.fastq.gz -2 ${filename2}.fastq.gz -p $num_threads --validateMappings -o $filename3.quant
          mv $filename3.quant $bam_out/Salmon_counts
        fi
      fi
    fi
}

paired_fastq()
{
    filename=$(basename "$f" ".fastq")
    filename2=${filename/_R1/_R2}
    filename3=$(echo $filename | sed 's/_R1//')
    if [ "$referencegenome" != 0 ] || [ "$index_folder" != 0 ]; then
      if [ "$transcriptome" == 0 ]; then
        if [ ! -d "$bam_out" ]; then
              mkdir $bam_out
        fi
        if [ "$fastqc" != 0 ]; then
          DIRECTORY=$bam_out/Fastqc_out
          if [ ! -d "$DIRECTORY" ]; then
            mkdir $DIRECTORY
          fi
          fastqc ${filename}.fastq ${filename2}.fastq
          mkdir "$filename3".fastqc_out
          mv *.zip *.html "$filename3".fastqc_out
          mv "$filename3".fastqc_out $DIRECTORY
        fi
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename3 "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Mapping step"
        echo "#######################"
        if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
        elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename3.sam --met-file metrics.txt 2>> mapped.txt
        fi    
        echo "samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam"
        samtools view -@ $num_threads -Sbo $filename3.bam $filename3.sam
        duplicates_paired
        stringtie_non_SRA
        if [ "$hisat" != 0 ]; then
          featurecounts
        fi
        rm $filename3.sam
      else
        if [ "$transcriptome" != 0 ]; then
          mkdir -p $bam_out/Salmon_counts
          if [ "$fastqc" != 0 ]; then
            DIRECTORY=$bam_out/Fastqc_out
            if [ ! -d "$DIRECTORY" ]; then
              mkdir $DIRECTORY
            fi
            fastqc ${filename}.fastq ${filename2}.fastq
            mkdir "$filename3".fastqc_out
            mv *.zip *.html "$filename3".fastqc_out
            mv "$filename3".fastqc_out $DIRECTORY
          fi
          salmon quant -i salmon_index -l A -1 ${filename}.fastq -2 ${filename2}.fastq -p $num_threads --validateMappings -o $filename3.quant
          mv $filename3.quant $bam_out/Salmon_counts
        fi
      fi
    fi
}

single_end()
{
    extension=$(echo "$f" | sed -r 's/.*(fq|fq.gz|fastq|fastq.gz)$/\1/')
    filename=$(basename "$f" ".$extension")
    if [ "$referencegenome" != 0 ] || [ "$index_folder" != 0 ]; then
      if [ "$transcriptome" == 0 ]; then
        if [ ! -d "$bam_out" ]; then
              mkdir $bam_out
        fi
        echo "------------------------------------------" >> mapped.txt
        echo "### Mapping percentages of" $filename "###" >> mapped.txt
        echo "------------------------------------------" >> mapped.txt
        echo "#######################"
        echo "Running Mapping step"
        echo "#######################"
        if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
          if [ "$lib_type" == "US" ]; then
            echo "hisat2 -x $fbname -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
          else
            echo "hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
            hisat2 -x $fbname --rna-strandness $lib_type -U $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
          fi
        elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
            echo "bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt"
            bowtie2 -x ref_genome -U $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $filename.sam --met-file metrics.txt 2>> mapped.txt
        fi
        echo "samtools view -@ $num_threads -Sbo $filename.bam $filename.sam"
        samtools view -@ $num_threads -Sbo $filename.bam $filename.sam
        echo "sambamba sort --tmpdir=temp -t $num_threads -o $filename.sorted.bam $filename.bam"
        sambamba sort --tmpdir=temp -t $num_threads -o $filename.sorted.bam $filename.bam
        rm -r temp $filename.bam
        if [ "$duplicates" != 0 ]; then
          rmdup=$(basename $filename.sorted.bam ".sorted.bam")
          picard MarkDuplicates I=$filename.sorted.bam O=$rmdup."sorted.rmdup.bam" ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
          mv $filename.sorted.bam $filename.sorted.bam.bai $rmdup".sorted.rmdup.bam" "$bam_out"
        else
          mv $filename.sorted.bam $filename.sorted.bam.bai "$bam_out"
        fi
        stringtie_non_SRA_single
        if [ "$hisat" != 0 ]; then
          featurecounts
        fi
        rm $filename.sam
      else
        if [ "$transcriptome" != 0 ]; then
           mkdir -p $bam_out/Salmon_counts
      	   salmon quant -i salmon_index -l A -r $f -p $num_threads --validateMappings -o $f.quant
           mv $f.quant $bam_out/Salmon_counts
        fi
      fi
    fi
}

# ############################################################################################################################################################################################################################
# # Reference genome/Index ###
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
elif [ "$transcriptome" != 0 ]; then
  echo "####################################"
  echo "Starting transcriptome indexing"
  echo "####################################"
  echo "salmon index -i salmon_index -t $transcriptome"
  salmon index -i salmon_index -t "$transcriptome"
fi

# ############################################################################################################################################################################################################################
# # Mapping
# ############################################################################################################################################################################################################################

# Paired end reads

if [ ! -z "$left_reads" ] && [ ! -z "$right_reads" ]; then
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
    if [ "$transcriptome" != 0 ]; then
      exit
    fi
    success_message
  
# Single end reads

elif [ ! -z "$single_reads" ]; then
	numb=$(ls "${single_reads[@]}" | wc -l)
	for f in "${single_reads[@]}"; do
    if [ "$transcriptome" != 0 ]; then
      if [ ! -d "$bam_out" ]; then
            mkdir $bam_out
      fi
    fi
		if [ "$fastqc" != 0 ]; then
      DIRECTORY=$bam_out/Fastqc_out
      if [ ! -d "$DIRECTORY" ]; then
        mkdir -p $DIRECTORY
      fi
  		fastqc $f
  		mkdir "$f".fastqc_out
  		mv *.zip *.html "$f".fastqc_out
  		mv "$f".fastqc_out $DIRECTORY
    fi
    single_end
  done
    if [ "$transcriptome" == 0 ]; then
    	success_message
    elif [ "$transcriptome" != 0 ]; then
		  exit
	  fi
  
# SRA

elif [ ! -z $sra_id ]; then
  if [[ -f $sra_id ]]; then
  	if [ "$referencegenome" != 0 ] || [ "$index_folder" != 0 ]; then
      if [ "$transcriptome" == 0 ]; then
  	      mkdir "$bam_out" 
    	    while read f; do
  	      echo "------------------------------------------" >> mapped.txt
  	      echo "### Mapping percentages of" $f "###" >> mapped.txt
  	      echo "------------------------------------------" >> mapped.txt
  	      echo "#######################"
  	      echo "Running Mapping step"
  	      echo "#######################"
  	      if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then #(-Q -h)
  	        if [ "$lib_type" == "US" ]; then 
  	          echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt
  	        else
  	          echo "hisat2 -x $fbname --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt
  	        fi
  	      elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then #(-Q -b)
  	          echo "bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
  	          bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $f.sam --met-file metrics.txt 2>> mapped.txt 
  	      elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default #(-h)
  	        if [ "$lib_type" == "US" ]; then
  	          echo "hisat2 -x $fbname --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt
  	        else
  	          echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt
  	        fi
  	      elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then #(-b)
  	          echo "bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt"
  	          bowtie2 -x ref_genome --sra-acc $f -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $f.sam --met-file metrics.txt 2>> mapped.txt
  	      fi
  	      samtools view -@ $num_threads -Sbo $f.bam $f.sam
  	      rm $f.sam
  	      echo "sambamba sort --tmpdir=temp -t $num_threads -o $f.sorted.bam $f.bam"
  	      sambamba sort --tmpdir=temp -t $num_threads -o $f.sorted.bam $f.bam
  	      if [ "$duplicates" != 0 ]; then
  	        rmdup=$(basename $f.sorted.bam ".sorted.bam")
  	        picard MarkDuplicates I=$f.sorted.bam O=$rmdup."sorted.rmdup.bam" ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
  	        mv $rmdup."sorted.rmdup.bam" "$bam_out"
  	      fi
  	      stringtie_SRA_multi
  	      if [ "$fastqc" != 0 ]; then
  	        if [ "$seq_type" == "SE" ]; then
  	          picard SamToFastq I="$f".bam FASTQ="$f".fastq
  	          fastqc "$f".fastq
  	          mkdir "$f".fastqc_out
  	          rm "$f".fastq
  	          mv *.zip *.html "$f".fastqc_out
  	        elif [ "$seq_type" == "PE" ]; then
  	          picard SamToFastq I="$f".bam FASTQ="$f"."R1".fastq SECOND_END_FASTQ="$f"."R2".fastq
  	          fastqc "$f"."R1".fastq "$f"."R2".fastq
  	          mkdir "$f".fastqc_out
  	          rm "$f".*.fastq
  	          mv *.zip *.html "$f".fastqc_out
  	        fi
  	      fi
  	    done < "$sra_id"
  	      if [ "$fastqc" != 0 ]; then
            DIRECTORY=$bam_out/Fastqc_out
            if [ ! -d "$DIRECTORY" ]; then
              mkdir -p $DIRECTORY
            fi
  	        mv *.fastqc_out $DIRECTORY
  	      fi
  	      if [ "$hisat" != 0 ]; then
  	        featurecounts
  	      fi
  	      rm -r temp *.bam
  	      success_message
	    else
    		while read f; do
    	      echo "#######################"
    	      echo "Running Salmon step"
    	      echo "#######################"
    			if [ "$transcriptome" != 0 ]; then
    				if [ "$seq_type" == "SE" ]; then
    					fastq-dump -A $f --gzip
      				salmon quant -i salmon_index -l A -r $f."fastq.gz" -p $num_threads --validateMappings -o $f.quant
              mkdir -p $bam_out/Salmon_counts
              mv $f.quant $bam_out/Salmon_counts
      				if [ "$fastqc" != 0 ]; then
      					mkdir $f.fastqc
      					fastqc $f."fastq.gz" -o $f.fastqc
      					rm $f."fastq.gz"
      				else
      					rm $f."fastq.gz"
      				fi
      			elif [ "$seq_type" == "PE" ]; then
      				fastq-dump -A $f --gzip --split-files
      				salmon quant -i salmon_index -l A -1 $f"_1".fastq.gz -2 $f"_2".fastq.gz -p $num_threads --validateMappings -o $f.quant
              mkdir -p $bam_out/Salmon_counts
              mv $f.quant $bam_out/Salmon_counts
      				if [ "$fastqc" != 0 ]; then
      					mkdir $f"_1".fastqc $f"_2".fastqc
      					fastqc $f"_1".fastq.gz -o $f"_1".fastqc
      					fastqc $f"_2".fastq.gz -o $f"_2".fastqc
      					rm $f"_1".fastq.gz $f"_2".fastq.gz
      				else
      					rm $f"_1".fastq.gz $f"_2".fastq.gz
      				fi
          	fi
          fi
    	 done < "$sra_id"
      if [ "$fastqc" != 0 ]; then
        DIRECTORY=$bam_out/Fastqc_out
        if [ ! -d "$DIRECTORY" ]; then
          mkdir -p $DIRECTORY
        fi
        mv *.fastqc $DIRECTORY
      fi
      echo "##############################"
      echo "Pipeline executed successfully"
      echo "##############################"
    fi  
    fi
  else
  	if [ "$referencegenome" != 0 ] || [ "$index_folder" != 0 ]; then
      if [ "$transcriptome" == 0 ]; then
      	mkdir "$bam_out"
      	echo "------------------------------------------" >> mapped.txt
      	echo "### Mapping percentages of" $sra_id "###" >> mapped.txt
      	echo "------------------------------------------" >> mapped.txt
      	echo "#######################"
      	echo "Running Mapping step"
      	echo "#######################"
          if [ "$quality_64" != 0 ] && [ "$hisat" != 0 ]; then
  	        if [ "$lib_type" == "US" ]; then 
  	          echo "hisat2 -x $fbname --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
  	        else
  	          echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
  	        fi
  	    elif [ "$quality_64" != 0 ] &&[ "$bowtie" != 0 ]; then
  	        echo "bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
  	        bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred64 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
  	    elif [ "$quality_64" == 0 ] && [ "$hisat" != 0 ]; then # phred 33 default
  	        if [ "$lib_type" == "US" ]; then 
  	          echo "hisat2 -x $fbname --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
  	        else
  	          echo "hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
  	          hisat2 -x $fbname --rna-strandness $lib_type --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --min-intronlen $min_intl --max-intronlen $max_intl --dta --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
  	        fi
  	    elif [ "$quality_64" == 0 ] &&[ "$bowtie" != 0 ]; then
  	        echo "bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt"
  	        bowtie2 -x ref_genome --sra-acc $sra_id -p $num_threads -5 $five_trim -3 $three_trim --phred33 -S $sra_id.sam --met-file metrics.txt 2>> mapped.txt
  	    fi
  	    samtools view -@ $num_threads -Sbo $sra_id.bam $sra_id.sam
  	    rm $sra_id.sam
  	    echo "sambamba sort --tmpdir=temp -t $num_threads -o $sra_id.sorted.bam $sra_id.bam"
  	    sambamba sort --tmpdir=temp -t $num_threads -o $sra_id.sorted.bam $sra_id.bam
  	    if [ "$duplicates" != 0 ]; then
  	      rmdup=$(basename $sra_id.sorted.bam ".sorted.bam")
  	      picard MarkDuplicates I=$sra_id.sorted.bam O=$rmdup."sorted.rmdup.bam" ASSUME_SORTED=TRUE METRICS_FILE=/dev/null VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
  	      mv $rmdup."sorted.rmdup.bam" "$bam_out"
  	    fi
  	    stringtie_SRA_single
  	    if [ "$fastqc" != 0 ]; then
  	      if [ "$seq_type" == "SE" ]; then
  	        picard SamToFastq I=$sra_id.bam FASTQ="$sra_id".fastq
  	        fastqc "$sra_id".fastq
  	        mkdir fastqc_out
  	        rm "$sra_id".fastq
  	        mv *.zip *.html fastqc_out
  	        mv fastqc_out "$bam_out/Fastqc_out"
  	      elif [ "$seq_type" == "PE" ]; then
  	        picard SamToFastq I=$sra_id.bam FASTQ="$sra_id"."R1".fastq SECOND_END_FASTQ="$sra_id"."R2".fastq
  	        echo "fastqc "$sra_id"."R1".fastq "$sra_id"."R2".fastq"
  	        fastqc "$sra_id"."R1".fastq "$sra_id"."R2".fastq
  	        mkdir fastqc_out
  	        rm "$sra_id".*.fastq
  	        mv *.zip *.html fastqc_out
  	        mv fastqc_out "$bam_out/Fastqc_out"
  	      fi
  	    fi
  	    if [ "$hisat" != 0 ]; then
  	      featurecounts
  	    fi
  	    rm -r temp $sra_id.bam
  	    success_message
  	  else
		    if [ "$transcriptome" != 0 ]; then
           mkdir $bam_out
			     if [ "$seq_type" == "SE" ]; then
				     fastq-dump -A $sra_id --gzip
      	     salmon quant -i salmon_index -l A -r $sra_id."fastq.gz" -p $num_threads --validateMappings -o $sra_id.quant
             mkdir -p $bam_out/Salmon_counts
             mv $sra_id.quant $bam_out/Salmon_counts
  			       if [ "$fastqc" != 0 ]; then
                  DIRECTORY=$bam_out/Fastqc_out
                  if [ ! -d "$DIRECTORY" ]; then
                    mkdir -p $DIRECTORY
                  fi
  				        mkdir $sra_id.fastqc
  				        fastqc $sra_id."fastq.gz" -o $sra_id.fastqc
  				        rm $sra_id."fastq.gz"
                  mv $sra_id.fastqc $bam_out/Fastqc_out
                  echo "##############################"
                  echo "Pipeline executed successfully"
                  echo "##############################"
  			       else
  				        rm $sra_id."fastq.gz"
                  echo "##############################"
                  echo "Pipeline executed successfully"
                  echo "##############################"
  			       fi
  		     elif [ "$seq_type" == "PE" ]; then
  			     fastq-dump -A $sra_id --gzip --split-files
  			     salmon quant -i salmon_index -l A -1 $sra_id"_1".fastq.gz -2 $sra_id"_2".fastq.gz -p $num_threads --validateMappings -o $sra_id.quant
             mkdir -p $bam_out/Salmon_counts
             mv $sra_id.quant $bam_out/Salmon_counts
  			     if [ "$fastqc" != 0 ]; then
                DIRECTORY=$bam_out/Fastqc_out
                if [ ! -d "$DIRECTORY" ]; then
                  mkdir -p $DIRECTORY
                fi
      				  mkdir $sra_id"_1".fastqc $sra_id"_2".fastqc
      				  fastqc $sra_id"_1".fastq.gz -o $sra_id"_1".fastqc
      				  fastqc $sra_id"_2".fastq.gz -o $sra_id"_2".fastqc
      				  rm $sra_id"_1".fastq.gz $sra_id"_2".fastq.gz
                mv $sra_id"_1".fastqc $sra_id"_2".fastqc $bam_out/Fastqc_out
                echo "##############################"
                echo "Pipeline executed successfully"
                echo "##############################"
  			     else
  				      rm $sra_id"_1".fastq.gz $sra_id"_2".fastq.gz
                echo "##############################"
                echo "Pipeline executed successfully"
                echo "##############################"
  			     fi
      	   fi
        fi
  	  fi
    fi
  fi
fi
##### End ########
