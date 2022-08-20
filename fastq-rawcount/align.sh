#!/bin/sh

#  align.sh
#  
#
#  Created by Ken Nakatsu on 7/27/22.
#  

#With eight cores, this pipeline runs at a speed of 45 seconds per sample. ~4mil rd per sample

#this script will do:
# 1) run a quality check
# 2) remove adapters based upon a given adapter seq
# 3) run a post-trim quality check
# 4) run alignment with bowtie 2 against hg38
# 5) run annotation with four different files
# 5.1) all, all_smrna, mirna, pirna
#start the script in a dir with two files, annotation and unaligned files with ref. genome
#ensure that you have conda installed with all packages
eval "$(conda shell.bash hook)"
conda activate smrnaseq_ddocent_env

mv $unalign_dir unaligned

mkdir output
#all output files will be stored here
mkdir output/fastqc_analysis
mkdir output/fastqc_analysis/unprocessed
#fastq
mkdir output/fastqc_analysis/docs
#zips and logs
mkdir output/fastqc_analysis/results
#html files

#run fastqc with 8 threads. move to respective directory
fastqc unaligned/*.fastq.gz -t $threads -o output/fastqc_analysis
mv output/fastqc_analysis/*.html output/fastqc_analysis/results
mv output/fastqc_analysis/*.zip output/fastqc_analysis/docs

#user can input a new adapter or add another one with --adapter2 '$user_adapter2'
#this will remove adapters and move it to its respective dir

cd unaligned

#var%%.x will remove everything after . i.e. y.x -> y

for f in *.fastq.gz ; do
    file="$f"
    sampname="${file%%.*}"
    AdapterRemoval --file1 "${file}" --adapter1 $user_adapter1 --basename "trimmed"_"${sampname}" --gzip --threads $threads --trimqualities --minlength $minlen
done

cd ..

mkdir output/adaptertrim
mkdir output/adaptertrim/discarded_ar
mkdir output/adaptertrim/settings

mv unaligned/*.discarded.gz output/adaptertrim/discarded_ar
mv unaligned/*.settings output/adaptertrim/settings
mv unaligned/*.truncated.gz output/adaptertrim

#change into fastq file
#note: ${var//x/y} will take the variable, remove x from it, and replace it with y
for f in output/adaptertrim/*.gz ; do
    mv $f ${f//truncated/fastq};
done
#passes ~800,000-900,000 reads per second

#run post trimming QC
mkdir output/fastqc_post_trim
fastqc output/adaptertrim/*.fastq.gz -t $threads -o output/fastqc_post_trim

#now lets align and move those aligned files to a new folder

mkdir output/aligned

for f in output/adaptertrim/*.fastq.gz; do
    file="$f"
    sampname="${file%%.*}"
    bowtie2 -x unaligned/grch38_1kgmaj -U $f -S "${sampname}".sam -p $threads
done
#40% alignment for media
#90-98% for cell

for f in output/adaptertrim/*.sam; do
    mv $f output/aligned
done

#run featureCounts for each one of the annotation files

if (($reg==1))
then
    featureCounts -T $threads -a annotation/hsa_all.gff3 -F 'GTF' -g 'gene_name' -o "hsa_all.tsv" output/aligned/*.sam
fi
#~30%

if (($sm_all==1))
then
    featureCounts -T $threads -a annotation/hsa_all_smrna.gff -F 'GTF' -t 'lnc_RNA' -g 'gene_id'  -o "hsa_all_smrna.tsv" output/aligned/*.sam
fi
#30-40%

if (($mirna==1))
then
    featureCounts -T $threads -a annotation/hsa_mirna.gff3 -F 'GTF' -t 'miRNA' -g 'Name' -o "hsa_mirna.tsv" output/aligned/*.sam
fi
#7-50%

if (($pirna==1))
then
    featureCounts -T $threads -a annotation/hsa_pirna.gtf -F 'GTF' -t 'piRNA' -g 'piRNA_code' -o "hsa_pirna.tsv" output/aligned/*.sam
fi
#1-12%

mkdir output/counts
mv *.tsv output/counts
for f in *.summary ; do
    mv $f ${f//.tsv.summary/.summary.tsv};
done
mv *.tsv output/counts
mkdir output/counts/summary
mv output/counts/*.summary.tsv output/counts/summary

mkdir ../rawcount-analysis/counts
mv output/counts/*.tsv ../rawcount-analysis/counts

touch output/counts/note.txt
echo 'raw counts are now in the rawcount-analysis/counts folder' > output/counts/note.txt
#With eight cores, this pipeline runs at a speed of 45 seconds per sample, 18 minutes for 23 samples. 
