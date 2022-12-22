# EZanalyze (wip)
## Ken Nakatsu in Collaboration with the Deng Lab at John A. Burns School of Medicine

## Usage
1) Run scripts in installer, depending upon your system and what is already installed.
2) Run the scripts in the fastq-rawcount folder, ensuring to modify config.sh to your desire.
3) Once config.sh has been run, please run align.sh to generate aligned files. 
4) Then run the deseq2_analysis.R to obtain DEGs from your files. Configuration for the analysis part is still a WIP. However, the user, if familiar with R, should be able to modify the code to fit their samples.
5) Run organize.sh to store outputs nicely in folders.

## To do
- Commit new alignment script which properly deduplicates with UMIs, allowing for user input for the format of UMI. 
- Have alignement, count, and analyses run in one command. 
- Create additionally functionality in analysis of data. 
- Introduce functionality that produces a pdf output that nicely summarizes experiments with preliminary conclusions.
- To expand to other species, create an additional component that analyzes the input annotation files. (i.e. how much they overlap) as, often times, small ncRNA annotation databases overlap with others leading to misleading/ambiguous results. 


## Citations

## Bioconda
Grüning, Björn, Ryan Dale, Andreas Sjödin, Brad A. Chapman, Jillian Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda Team, and Johannes Köster. 2018. Bioconda: Sustainable and Comprehensive Software Distribution for the Life Sciences. Nature Methods, 2018 doi:10.1038/s41592-018-0046-7.

## FastQC
Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

## AdapterRemoval
Lindgreen (2012): AdapterRemoval: Easy Cleaning of Next Generation Sequencing Reads, BMC Research Notes, 5:337 http://www.biomedcentral.com/1756-0500/5/337/

## Bowtie2
Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

## FeatureCounts
Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014 Apr 1;30(7):923-30. doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PMID: 24227677.
## piRNA Database
 https://www.pirnadb.org/about/informations/pirna
## all smRNA Database
https://dashr2.lisanwanglab.org/about.php
##miRBase
miRBase: from microRNA sequences to function.
Kozomara A, Birgaoanu M, Griffiths-Jones S.
Nucleic Acids Res 2019 47:D155-D162

miRBase: annotating high confidence microRNAs using deep sequencing data.
Kozomara A, Griffiths-Jones S.
Nucleic Acids Res 2014 42:D68-D73

miRBase: integrating microRNA annotation and deep-sequencing data.
Kozomara A, Griffiths-Jones S.
Nucleic Acids Res 2011 39:D152-D157

miRBase: tools for microRNA genomics.
Griffiths-Jones S, Saini HK, van Dongen S, Enright AJ.
Nucleic Acids Res 2008 36:D154-D158

miRBase: microRNA sequences, targets and gene nomenclature.
Griffiths-Jones S, Grocock RJ, van Dongen S, Bateman A, Enright AJ.
Nucleic Acids Res 2006 34:D140-D144

The microRNA Registry.
Griffiths-Jones S.
Nucleic Acids Res 2004 32:D109-D111

## Regular hg38 annotation
##### https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/

