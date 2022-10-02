# Overview
The MNaseSeq github repository is stored locally, and will be used for project deployment. Multiple projects can be deployed from this one point simultaneously, without concern.

## 1. Getting Started
### 1.1 Introduction¶
The MNaseSeq Pipeline beings with raw FASTQ files and performs adaptor trimming, assembly, and alignment. Bed files are created, and depending on user input, selected regions of interst may be used. Fragment centers (DYAD's) are then determined, and histograms of occurences are created. QC reports are also generated with each project.

The following are sub-commands used within MNaseSeq:
- init: initalize the pipeline
- dryrun: predict the binding of peptides to any MHC molecule
- run: execute the pipeline on the Biowulf HPC
- runlocal: execute a local, interactive, session
- unlock: unlock directory
- reset: delete a workdir, and re-initialize

## 1.2 Setup Dependencies
MNaseSeq has several dependencies listed below. These dependencies can be installed by a sysadmin. All dependencies will be automatically loaded if running from Biowulf.

- bedtools: "bedtools/2.30.0"
- bowtie2: "bowtie/2-2.4.2"
- cutadapt: "cutadapt/1.18"
- pear: "pear/0.9.11"
- python: "python/3.7"
- R: "R/4.0.3"
- samtools: "samtools/1.11"

## 1.3 Login to the cluster
MNaseSeq has been exclusively tested on Biowulf HPC. Login to the cluster's head node and move into the pipeline location.
```
# ssh into cluster's head node
ssh -Y $USER@biowulf.nih.gov
```

## 1.4 Load an interactive session
An interactive session should be started before performing any of the pipeline sub-commands, even if the pipeline is to be executed on the cluster.
```
# Grab an interactive node
srun -N 1 -n 1 --time=12:00:00 -p interactive --mem=8gb  --cpus-per-task=4 --pty bash
```