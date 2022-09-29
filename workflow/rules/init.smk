#########################################################
# IMPORT PYTHON LIBRARIES HERE
#########################################################
import sys
import os
import pandas as pd
import yaml
import subprocess
# import glob
# import shutil
#########################################################


#########################################################
# FILE-ACTION FUNCTIONS 
#########################################################
def check_existence(filename):
  if not os.path.exists(filename):
    exit("# File: %s does not exists!"%(filename))

def check_readaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.R_OK):
    exit("# File: %s exists, but cannot be read!"%(filename))

def check_writeaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.W_OK):
    exit("# File: %s exists, but cannot be read!"%(filename))

def get_file_size(filename):
    filename=filename.strip()
    if check_readaccess(filename):
        return os.stat(filename).st_size
#########################################################

#########################################################
# DEFINE CONFIG FILE AND READ IT
#########################################################
CONFIGFILE = str(workflow.overwrite_configfiles[0])

# set memory limit 
# used for sambamba sort, etc
# MEMORYG="100G"

# read in various dirs from config file
WORKDIR=config['workdir']
RESULTSDIR=join(WORKDIR,"results")

# get scripts folder
SCRIPTSDIR = join(WORKDIR,"scripts")
check_existence(SCRIPTSDIR)

# get resources folder
try:
    RESOURCESDIR = config["resourcesdir"]
except KeyError:
    RESOURCESDIR = join(WORKDIR,"resources")
check_existence(RESOURCESDIR)

# create symlink FQ dir
if not os.path.exists(join(WORKDIR,"fastqs")):
    os.mkdir(join(WORKDIR,"fastqs"))
if not os.path.exists(RESULTSDIR):
    os.mkdir(RESULTSDIR)

# check read access to required files
for f in ["samplemanifest"]:
    check_readaccess(config[f])
if config["run_select_bed"]=="Y":
    check_readaccess(config["intervals_of_interest"])
#########################################################

#########################################################
# CHECK MANIFESTS
#########################################################
# run script to check naming / fastq / metadata information is valid against requirements
python_script = join(SCRIPTSDIR,"check_manifest.py")
python_cmd= "python " + python_script + " " + RESULTSDIR + "/ " + config["samplemanifest"] + " " + config["contrastmanifest"]
subprocess.call(python_cmd, shell=True)
check_existence(join(RESULTSDIR,"manifest_qc_pass.txt"))
#########################################################

#########################################################
# CREATE SAMPLE DATAFRAME
#########################################################
# each line in the samplemanifest is a sample
SAMPLESDF = pd.read_csv(config["samplemanifest"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)

print("# Checking Sample Manifest...")
print("## \tTotal Samples in manifest : "+str(len(SAMPLES)))
print("# Checking read access to raw fastqs...")

# only PE samples are supported by the pipeline
for sampleid in SAMPLES:
    # set FASTQ files
    R1file=SAMPLESDF["path_to_R1_fastq"][sampleid]
    R2file=SAMPLESDF["path_to_R2_fastq"][sampleid]

    # check R1 access
    check_readaccess(R1file)
    R1filenewname=join(WORKDIR,"fastqs",sampleid+".R1.fastq.gz")
    if not os.path.exists(R1filenewname):
        os.symlink(R1file,R1filenewname)
        SAMPLESDF.loc[[sampleid],"R1"]=R1filenewname
    
    # check $2 access
    check_readaccess(R2file)
    R2filenewname=join(WORKDIR,"fastqs",sampleid+".R2.fastq.gz")
    if not os.path.exists(R2filenewname):
        os.symlink(R2file,R2filenewname)
        SAMPLESDF.loc[[sampleid],"R2"]=R2filenewname

print("## \tRead access to all raw fastqs is confirmed!")

#########################################################

#########################################################
# READ IN TOOLS REQUIRED BY PIPELINE
# THESE INCLUDE LIST OF BIOWULF MODULES (AND THEIR VERSIONS)
# MAY BE EMPTY IF ALL TOOLS ARE DOCKERIZED
#########################################################
## Load tools from YAML file
try:
    TOOLSYAML = config["tools"]
except KeyError:
    TOOLSYAML = join(WORKDIR,"resources","tools.yaml")
check_readaccess(TOOLSYAML)
with open(TOOLSYAML) as f:
    TOOLS = yaml.safe_load(f)
#########################################################

#########################################################
# READ CLUSTER PER-RULE REQUIREMENTS
#########################################################
## Load cluster.json
try:
    CLUSTERYAML = config["clusteryaml"]
except KeyError:
    CLUSTERYAML = join(WORKDIR,"resources","cluster.yaml")
check_readaccess(CLUSTERYAML)
with open(CLUSTERYAML) as yaml_file:
    CLUSTER = yaml.load(yaml_file, Loader=yaml.FullLoader)

## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")
#########################################################

#########################################################
# QC parameters
#########################################################
QCDIR=join(RESULTSDIR,"QC")

#FASTQ_SCREEN_CONFIG=config["fastqscreen_config"]
#check_readaccess(FASTQ_SCREEN_CONFIG)
#print("# FQscreen config  :",FASTQ_SCREEN_CONFIG)
#########################################################

#########################################################
# SET OTHER PIPELINE GLOBAL VARIABLES
#########################################################
print("#"*66)
print("#"*66)
print("# Pipeline Parameters:")
print("# Working dir :",WORKDIR)
print("# Results dir :",RESULTSDIR)
print("# Scripts dir :",SCRIPTSDIR)
print("# Resources dir :",RESOURCESDIR)
print("# Cluster YAML :",CLUSTERYAML)

#GENOME=config["genome"]
#INDEXDIR=config[GENOME]["indexdir"]
#print("# Bowtie index dir:",INDEXDIR)

#GENOMEFILE=join(INDEXDIR,GENOME+".genome") # genome file is required by macs2 peak calling
#check_readaccess(GENOMEFILE)
#print("# Genome :",GENOME)
#print("# .genome :",GENOMEFILE)

#GENOMEFA=join(INDEXDIR,GENOME+".fa") # genome file is required by motif enrichment rule
#check_readaccess(GENOMEFA)
#print("# Genome fasta:",GENOMEFA)
#########################################################
