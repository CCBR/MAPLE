# 2. Preparing Files
The pipeline is controlled through editing configuration and manifest files. Defaults are found in the /WORKDIR/ after initialization.

## 2.1 Configs
The configuration files control parameters and software of the pipeline. These files are listed below:

- resources/cluster.yaml
- resources/tools.yaml
- config.yaml

### 2.1.1 Cluster YAML
The cluster configuration file dictates the resouces to be used during submission to Biowulf HPC. There are two differnt ways to control these parameters - first, to control the default settings, and second, to create or edit individual rules. These parameters should be edited with caution, after significant testing.

### 2.1.2 Tools YAML
The tools configuration file dictates the version of each tool that is being used. Updating the versions may break specific rules if versions are not backwards compatible with the defaults listed.

### 2.1.3 FASTQ Screen Config
The FASTQ screen configuration files dictates the parameters used for FASTQ Screen.

### 2.1.4 MultiQC Config
The MultiQC screen configuration files control the parameters used for MultiQC.

### 2.1.5 Config YAML
There are several groups of parameters that are editable for the user to control the various aspects of the pipeline. These are :

- Folders and Paths
      - These parameters will include the input and ouput files of the pipeline, as well as list all manifest names.
- User parameters
      - These parameters will control the pipeline features. These include thresholds and whether to perform processes.

## 2.2 Preparing Manifests
There are two manifests which are required for the pipeline. These files describe information on the samples and desired contrasts. The paths of these files are defined in the config.yaml file. These files are:

- sampleManifest
- contrastManifest

### 2.2.1 Samples Manifest
This manifest will include information to sample level information. It includes the following column headers: sampleName type path_to_R1_fastq path_to_R2_fastq

- sampleName: the sampleID associated with the fasta file; which are unique. This may be a shorthand name, and will be used throughout the analysis.
- type: demographic information regarding the sample; example 'tumor'
- path_to_R1_fastq: the full path to the R1.fastq.gz file
- path_to_R1_fastq: the full path to the R2.fastq.gz file

An example sampleManifest file with multiplexing of one sample. Notice that the multiplexID test_1 is repeated, as Ro_Clip and Control_Clip are both found in the same fastq file, whereas test_2 is not multiplexed:

```
sampleName  type    path_to_R1_fastq                path_to_R2_fastq
Sample1     tumor   /path/to/sample1.R1.fastq.gz    /path/to/sample1.R2.fastq.gz
Sample2     tumor   /path/to/sample2.R1.fastq.gz    /path/to/sample2.R2.fastq.gz
Sample3     tumor   /path/to/sample3.R1.fastq.gz    /path/to/sample3.R2.fastq.gz
Sample4     tumor   /path/to/sample4.R1.fastq.gz    /path/to/sample4.R2.fastq.gz
```

### 2.2.2 Contrast Manifest
This manifest will include sample level information to compare samples. - contrast1: the sample name, identified in the samplesManifest [sampleName] column, of the sample to compare. example: 'Sample1' - contrast2: the sample name, identified in the samplesManifest [sampleName] column, of the sample to compare. example: 'Sample2'

Additional columns may be added, however, the same number of comparisons must be performed each time. For example, manifest 1 would be accepted, but manifest 2 would fail, as row 1 has three comparisons and row 2 has four comparisons.

Manifest example 1 (PASS)
```
contrast1   contrast2   contrast3
Sample1     Sample2     Sample3
Sample2     Sample3     Sample4
```

Manifest example 2 (FAIL)
```
contrast1   contrast2   contrast3   contrast4
Sample1     Sample2     Sample3
Sample2     Sample3     Sample4     Sample1
```