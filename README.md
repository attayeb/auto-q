# Auto-q
Qiime Analysis Automating Script.

This script is written to reduce the effort and time for Qiime analysis.

## Installation:
This script is designed to be installed in Qiime virtual machine. first copy the file to specific folder and then add it to your path.

After that you need to modify the configuration file (qiime.cfg) to update the folders of databases according to your settings.

## Steps of analysis:
You need to prepare Fastq files in one folder, not compressed:
```
usage: auto-q.py [-h] -i INPUT -o OUTPUT [-b BEGINWITH] [-t TRIM_THRESHOLD]
                 [-s STOP_AT] [-j JOINING_METHOD] [-d FASTQ_P]
                 [-q QC_THRESHOLD] [-c CONFIGFILE]
                 [--adapter ADAPTER_REFERENCE] [-a MAPPING_FILE]
                 [-p PARAMETER_FILE_NAME] [-n NUMBER_OF_CORES] [-f] [-m]
                 [-r RDB] [-e DEPTH] [--ml MINIMUM_LENGTH]

```


```
optional arguments:
  -h, --help            show this help message and exit
  -i INPUT              The folder where Fastq files are stored [required]
  -o OUTPUT             The folder of all results [required]
  -b BEGINWITH          begin with: [otu_picking], [diversity_analysis]
  -t TRIM_THRESHOLD     phred quality threshold for trimming [10]
  -s STOP_AT            stop at [chimera_removal\]
  -j JOINING_METHOD     choose the merging method (fastq-join) or (bbmerge)
  -d FASTQ_P            Percentage of maximum difference in fastq-join [16]
  -q QC_THRESHOLD       quality control phred threshold [19]
  -c CONFIGFILE         Configuration file name [qiime.cfg]
  --adapter ADAPTER_REFERENCE
                        Adapters reference file
  -a MAPPING_FILE       Mapping file name
  -p PARAMETER_FILE_NAME
                        The name of the parameter file [if not assigned is
                        automatically produced using configuration file
  -n NUMBER_OF_CORES    Number of cores to be used for the analysis [2]
  -f                    Using Unite database for fungal samples[False]
  -m                    Assign maxloose to be true for bbmerge [False]
  -r RDB                Reference data base [silva, greengenes]
  -e DEPTH              set the depth of diversity analyses [10000]
  --ml MINIMUM_LENGTH   Minimum length of reads kept after merging [380]

```

## Example:
```
$ auto-q.py -i /data/experiment1/fastqs/ -o /data/experiment1/results/ -t 12 -d 10 -r silva -n 10 -e 5000 -c /bin/auto-q/qiime.cfg 
```

## Results:
Output folder will has 7 subfolders:

| Folder name | content                                   |
|-------------|-------------------------------------------|
| others\     | log file, Mapping file, parameter file    |
| trimmed\    | fastq files after trimming                |
| merged\     | fastq files after merging pair reads      |
| qc\         | fasta files after quality step            | 
| chi\        | fastq files after chimera removed         | 
| otus\       | picked otus *standard Qiime output*       |
| div\        | diversity analyses results                |

