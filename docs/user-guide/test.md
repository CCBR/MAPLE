# 5. Pipeline Tutorial
Welcome to the MNaseSeq Pipeline Tutorial!

## 5.1 Getting Started
Review the information on the Getting Started for a complete overview the pipeline. The tutorial below will use test data available on NIH Biowulf HPC only. All example code will assume you are running v1.0 of the pipeline, from the shared [tobedetermined] storage directory, using test_1 data.

**A. Change working directory to the iCLIP repository**
```
# general format
cd /path/to/pipeline/[version number]

# example
cd /path/to/pipeline/v1.0
```

**B. Initialize Pipeline**
```
./run --runmode=init --workdir=/path/to/output/dir
```

## 5.2 Prepare the test set
**A. Two different test data sets are available, depending on the need.**
These include:
- test_1: Single sample comparison (two samples)
- test_2: Triple sample comparison (three samples)

**B. Pull the test data to your output directory**

NOTE: Test data is currently available for v1.0. Please contact samantha.sevilla@nih.gov to create other test data.

```
# general format
sh /data/CCBR_Pipeliner/ccbr1214/test/run_test.sh \
    -t test_number \
    -v version_id \
    -o /path/to/output/dir

# example running test_1, v1.0:
sh /data/CCBR_Pipeliner/ccbr1214/test/run_test.sh \
    -t test_1 \
    -v v1.0 \
    -o /path/to/output/dir 
```

## 5.3 Complete dry-run
A. Complete a dry-run and review output
```
./run --runmode=dryrun --workdir=/path/to/output/dir
```
Ensure that an expected output is displayed. An expected output for test_1 is as follows:
```
job              count    min threads    max threads
-------------  -------  -------------  -------------
alignment            2             48             48
all                  1              1              1
assembly             2             48             48
dyad_analysis        2             48             48
hist_frags           2             48             48
trim_adaptors        2             48             48
total               11              1             48
```

An expected output for test_2 is as follows:
```
```

## 5.4 Run the pipeline
Execute pipeline on the cluster
```
#submit to the cluster
./run --runmode=run --workdir=/path/to/output/dir
```

## 5.5 Review outputs
Review the expected outputs on the Output page. If there are errors, review and performing stesp described on the Troubleshooting page as needed.

##  5.6 Expert User Test Run
An example of a test workflow is shown below for test_1, v1.0 of the pipeline
```
source_dir="/home/sevillas2/git/ccbr1214";\
output_dir="/data/sevillas2/ccbr1214";\
test_id="test_1";\
version_id="v1.0";\
output_dir_full=${output_dir}/${version_id};\
if [[ -d $output_dir_full ]]; then rm -r $output_dir_full; fi ;\
cd $source_dir;\
./run --runmode=init --workdir=$output_dir_full;\
sh /data/CCBR_Pipeliner/ccbr1214/test/run_test.sh -t $test_id -v $version_id -o $output_dir_full;\
./run --runmode=dryrun --workdir=$output_dir_full;
./run --runmode=run --workdir=$output_dir_full;
```