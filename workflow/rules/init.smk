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
# CREATE CONTRAST DATAFRAME
#########################################################
if config["pipeline_phase"]=="third_pass":
    DACDF = pd.read_csv(config["contrastmanifest"],sep="\t",header=0)
    CONTRAST_FILES = DACDF['DAC_files'].tolist()

    # pull the sample id from the contrasts DF to create final output file
    clean_list=list()
    for contID in CONTRAST_FILES:
        check_readaccess(contID)
        save=contID.rsplit("/",1)[1]
        save=save.split(".",1)[0]
        clean_list.append(save)
    CONTRASTS_CLEAN_LIST='_AND_'.join(clean_list)
else:
    CONTRAST_FILES=config["samplemanifest"]
    CONTRASTS_CLEAN_LIST=""
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
# selected_bed file
#########################################################
# if user specifies to create a selected_bed
check_readaccess(config["master_bed_file"])

# set selection shorthand
bed_list=join(WORKDIR,"resources","bed_lists.csv",)
bed_df = pd.read_csv(bed_list)
selected_shorthand=bed_df['selected_shorthand'].tolist()
selected_bedfiles=bed_df['selected_bed'].tolist()

#########################################################

#########################################################
# fragment size
#########################################################
min_length=config["fragment_length_min"]
max_length=config["fragment_length_max"]

if not min_length.isnumeric() and max_length.isnumeric():
  raise Error("min_length and max_length must be integers")
else:
  if not max_length > min_length:
    raise Error("min_length must be larger than max_length")

#########################################################
# SET OTHER PIPELINE GLOBAL VARIABLES
#########################################################
print("# Pipeline Parameters:")
SPECIES=config["species"]
REF_SOURCE=config["reference_source"]
#########################################################

#########################################################
# SET INDEX VARIABLES
#########################################################
try:
    INDEXYAML = config["indexyaml"]
except KeyError:
    INDEXYAML = join(WORKDIR,"resources","index.yaml")
check_readaccess(INDEXYAML)
with open(INDEXYAML) as yaml_file:
    INDEX = yaml.load(yaml_file, Loader=yaml.FullLoader)

# check index access
INDEXDIR=config["index_dir"]
check_writeaccess(INDEXDIR)

# set flag to create index if it doesn't exist
GENOME_INDEX_FILE=join(INDEX[SPECIES][REF_SOURCE]["generated"] + ".1.bt2")
create_index_flag="N"
try:
    check_existence(GENOME_INDEX_FILE)
    print("## \tBowtie index source:",INDEX[SPECIES][REF_SOURCE]["generated"])
except:
    create_index_flag="Y"
    print("## \tBowtie index source will be created here:",INDEX[SPECIES][REF_SOURCE]["generated"])
INDEX_LOCATION = INDEX[SPECIES][REF_SOURCE]["generated"]
#########################################################

#########################################################
# Print values
#########################################################
print("## \tWorking dir :",WORKDIR)
print("## \tResults dir :",RESULTSDIR)
print("## \tScripts dir :",SCRIPTSDIR)
print("## \tResources dir :",RESOURCESDIR)
print("## \tCluster YAML :",CLUSTERYAML)
#########################################################
