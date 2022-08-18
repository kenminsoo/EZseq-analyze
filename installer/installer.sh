#!/bin/bash

#first install large packages
gdown_pip=$(pip show gdown)
gdown_state=$(echo ${test+1})
if [gdown_state=''];then
    pip install gdown
fi

#install anaconda3
cur_user=$(id -un)

bash Anaconda3-2022.05-MacOSX-x86_64.sh -b
source /Users/$cur_user/anaconda3/bin/activate

#conda update --all --yes

#install bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda create -n ezanalyze --yes

conda activate ezanalyze

#create separate environment for EZanalyze
#try to keep clean

conda install -c bioconda fastqc --yes
conda install -c bioconda adapterremoval --yes
conda install -c bioconda bowtie2 --yes
conda create -c bioconda -n featurecounts subread --yes
