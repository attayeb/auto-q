# QANAUS
Qiime 1 ANalysis AUtomating Script
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

