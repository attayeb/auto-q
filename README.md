# QANAUS
Qiime 1 ANalysis AUtomating Script.

This script is written to reduce the effort and time for Qiime 1 analysis.

## Installation:
This script is designed to be installed in Qiime 1 virtual machine. first copy the file to specific folder and then add it to your path.

After that you need to modify the configuration file (qiime.cfg) to update the folders of databases according to your settings.

## Steps of analysis:
You need to prepare Fastq files in one folder, not compressed:
```
usage: qiime_analysis.py [-h] -i INPUT -o OUTPUT [-s] [-f] [-t TRIM_THRESHOLD]
                         [-q QC_THRESHOLD] [-j] [-c CONFIGFILE]
                         [-p PARAMETER_FILE_NAME] [-n NUMBER_OF_CORES]
                         [-d FASTQ_P] [-m] [-e DEPTH]
qiime_analysis.py: error: argument -i is required
```


```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT              The folder where Fastq files are stored [required]
  -o OUTPUT             The folder of all results [required]
  -s                    Using SILVA database [False]
  -f                    Using Unite database for fungal samples[False]
  -t TRIM_THRESHOLD     phred quality threshold for trimming [12]
  -q QC_THRESHOLD       quality control phred threshold [19]
  -j                    use fastq-join for joining [False]
  -c CONFIGFILE         Configuration file name [qiime.cfg]
  -p PARAMETER_FILE_NAME
                        The name of the parameter file [if not assigned is
                        automatically produced using configuration file
  -n NUMBER_OF_CORES    Number of cores to be used for the analysis [2]
  -d FASTQ_P            Percentage of maximum difference in fastq-join [8]
  -m                    Assign maxloose to be true for bbmerge [False]
  -e DEPTH              set the depth of diversity analyses [10000]
```

## Example:
```buildoutcfg
$ qiime_analysis.py -i /data/experiment1/fastqs/ -o /data/experiment1/results/ -s -t 12 -p 10 -e 5000
```

## Results:
Output folder will has 5 folders:

| Folder name | content |
|--| -- |
| others | others|

|trimmed\| trimmed|

|merged\ |merged|
|qc\ |quality_step| 
|chi\| chimera_removed| 
|otus\| otus|
|div\ | diversity_analyses |