# Auto-q
Qiime Analysis Automating Script.

This script is written to reduce the effort and time for Qiime analysis.
It is designed to work on illumina pair-end reads FASTQ files.

## Installation:

1. Install usearch from <https://www.drive5.com/usearch/download.html> for a license issues you need to provide your email to install the free 32 bit version. QIIME works with version 6.1.544 please download it. 

Create bin/ folder in qiime home folder in the virtual machine /home/qiime/bin
```buildoutcfg
$ cd ~
$ mkdir ~/bin/
$ cd bin
$ mkdir usearch
```
copy usearch6.1.544_i86linux32 to /home/qiime/bin/usearch/ and rename the file name to usearch61

```buildoutcfg
$ mv usearch6.1.544_i86linux32 usearch61
```
To make usearch61 executable run the following command in usearch folder
```
$ chmod +x usearch61
``` 
[comment]: <> (this is comment)

2. Install bbduk tools from <https://sourceforge.net/projects/bbmap/> or you can run these commands in /home/qiime/bin/ folder: 

```buildoutcfg
$ wget https://sourceforge.net/projects/bbmap/files/BBMap_37.66.tar.gz
$ tar -zxvf BBMap_37.66.tar.gz && rm BBMap_37.66.tar.gz

```

This script is designed to be installed in Qiime virtual machine. first copy the file to specific folder and then add it to your path.

After that you need to modify the configuration file (qiime.cfg) to update the folders of databases according to your settings.

If you use windows and did not install Qiime virtual machine, please follow the instructions from Qiime web page: <http://qiime.org/install/virtual_box.html>

To run Auto-q, few softwares need to be installed in the Qiime virtual machine and need to be added to the execution path

3. Install Auto-q by running in home folder:
```buildoutcfg
$ git https://github.com/Attayeb/auto-q/ && rm -rf auto-q/.git 
$ cd auto-q
```
Edit .bashrc in your home directory and add the following line at the end:
```buildoutcfg
export PATH="/home/qiime/auto-q/:/home/qiime/bbtools/bin/:/home/qiime/usearch/:$PATH"
``` 

## Sequence files preparation:
### Fastq files: 
FASTQ files are named with the sample name and the sample number, which is a numeric assignment based on the order that the sample is listed in the sample sheet. Example:
                     
R1 &rarr;  SampleName_S1_L001_R1_001.fastq.gz 

R2 &rarr;  SampleName_S1_L001_R2_001.fastq.gz

keep a copy of the original compressed fastq files in a safe folder and use another copy after decompressing them. To decompress the fastq.gz file use this commnad inside the folder in terminal:
```
$ gunzip *.fastq.gz
``` 
Auto-q determines R1 and R2 using the name of the file, please do not modify them

## Steps of analysis:


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

## How to cite:
The paper of the script is under preparation now
