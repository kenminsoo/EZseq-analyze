#!/bin/bash

#have this file loaded in by a config option.

#usage: export var=value3
#1 = true
#0 = false
export unalign_dir=/Users/kenminsoo/Desktop/raw_data/*.fastq.gz
export multilane='no'
export threads=8
export user_adapter1="AACTGTAGGCACCATCAAT"
export user_adapter2=''
export minlen=19
export reg=1
export sm_all=1
export mirna=1
export pirna=1

bash align.sh
