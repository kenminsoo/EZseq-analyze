#!/bin/bash

#installer

#first let's get all the big files needed for installation and function
#python3 is needed
pipx install gdown

cur_user=$(id -un)

gdown https://drive.google.com/drive/folders/1uz8KoY63oWI_ayx_w7wjjlxkkkAinO_N -O /Users/$cur_user/Desktop/EZseq-analyze-main/installer  --folder
#change to find the actual location... this will do for now
#move alignment files to correct place.
mv *.bt2 ..
cd ..
mv *.bt2 fastq-rawcount/unaligned
cd installer

#install anaconda3

bash Anaconda3-2022.05-MacOSX-x86_64.sh -b
source /Users/$cur_user/anaconda3/bin/activate

#conda update --all --yes

#install bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n ezanalyze --yes

conda activate ezanalyze

#create separate environment for EZanalyze
#try to keep clean

conda install -c bioconda fastqc --yes
conda install -c bioconda adapterremoval --yes
conda install -c bioconda bowtie2 --yes
conda create -c bioconda -n featurecounts subread --yes
