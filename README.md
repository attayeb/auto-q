# Auto-q
Qiime Analysis Automating Script.

This script is written to reduce the effort and time for Qiime analysis.
It is designed to work on illumina pair-end reads FASTQ files.

## Installation:
1. This script is designed to be installed in Qiime virtual machine. To install QIIME, please check this link: <http://qiime.org/install/install.html>

2. Install usearch: follow the instructions in this link  <https://www.drive5.com/usearch/download.html> . QIIME works with version 6.1.544 32 bit. Please download it. 

* Create bin/ folder in qiime home folder in the virtual machine /home/qiime/bin
* Copy usearch6.1.544_i86linux32 to /home/qiime/bin/usearch/ and rename the file to usearch61
* Make usearch61 executable using this command
```
$ chmod +x usearch61
``` 

3. Install BBtools from <https://sourceforge.net/projects/bbmap/> or you can run these commands in /home/qiime/bin/ folder: 

```buildoutcfg
$ wget https://sourceforge.net/projects/bbmap/files/BBMap_37.66.tar.gz
$ tar -zxvf BBMap_37.66.tar.gz && rm BBMap_37.66.tar.gz

```

4. Install Auto-q by executing these commands in /home/qiime/bin/ :
```buildoutcfg
$ git clone https://github.com/Attayeb/auto-q/ && rm -rf auto-q/.git 
```
Edit .bashrc in your home directory and add the following line at the end:
```buildoutcfg
$ echo 'export PATH="/home/qiime/bin/auto-q/:/home/qiime/bin/bbtools/:/home/qiime/bin/usearch/:$PATH"' >> ~/.bashrc
```

5. If you want to use SILVA database you can download it from here <https://www.arb-silva.de/no_cache/download/archive/qiime/> use the latest one Silva_128_release.tgz, after downloading this file decompress it.
6. Modify qiime.cfg file to indicate the folders of your database. The default preinstalled greengenes folder is: /home/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/
 Please modify this file according to your settings.

## Sequence files preparation:
### Fastq files: 
FASTQ files are named with the sample name and the sample number, which is a numeric assignment based on the order that the sample is listed in the sample sheet. Example:
                     
R1 &rarr;  SampleName_S1_L001_R1_001.fastq.gz 

R2 &rarr;  SampleName_S1_L001_R2_001.fastq.gz

keep a copy of the original compressed fastq files in a safe folder and use another copy after 
decompressing them. To decompress the fastq.gz file use this commnad inside the folder in terminal:
```
$ gunzip *.fastq.gz
``` 
Auto-q determines R1 and R2 using the names of the files, please do not modify the file names.

## Steps of analysis:


```
usage: auto-q.py [-h] -i Input folder -o Output folder
                 [-t trim_phred_threshold] [-p fastq-join p]
                 [--adapter ADAPTER_REFERENCE] [-b starting step] [-s stop at]
                 [-j joining method] [-m] [-q quality control threshold]
                 [--continuation_reference newref_seq.fna]
                 [--continuation_otu_id C_OTU_ID] [-r Reference database]
                 [-c Configuration file name] [-a Mapping file name]
                 [--parameter_file_name PARAMETER_FILE_NAME]
                 [-n Number of jobs] [-e Sampling depth] [--ml Minimum length]

```


```
optional arguments:
  -h, --help            show this help message and exit
  -i Input folder       the input sequences filepath (fastq files) [REQUIRED]
  -o Output folder      the output directory [REQUIRED]
  -t trim_phred_threshold
                        phred quality threshold for trimming [default: 12]
  -p fastq-join p       fastq-join's percentage of mismatch [default: 16]
  --adapter ADAPTER_REFERENCE
                        Adapters reference file
  -b starting step      starting the analysis in the middle: (otu_picking),
                        (diversity_analysis)
  -s stop at            terminate the analysis at this step [choices:
                        (merging), (quality_control), (chimera_removal))
  -j joining method     choose the merging method (fastq-join) or (bbmerge)
                        [default: fastq-join]
  -m                    Assign maxloose to be true for bbmerge [default:
                        False]
  -q quality control threshold
                        quality control phred threshold [default: 19]
  --continuation_reference newref_seq.fna
                        reference sequence for continuation. If you want to
                        continue analysis using the reference data set from
                        previous analysis. you can find it in the last sample
                        otus folder new_refseqs.fna
  --continuation_otu_id C_OTU_ID
                        continuation reference new otus ids
  -r Reference database
                        silva, greengenes [default: silva]
  -c Configuration file name
                        Configuration file name [default: qiime.cfg]
  -a Mapping file name  Mapping file name
  --parameter_file_name PARAMETER_FILE_NAME
                        The name of the parameter file [if not assigned is
                        automatically produced using configuration file
  -n Number of jobs     Specify the number of jobs to start with [default: 2]
  -e Sampling depth     sampling depth for diversity analyses [default: 10000]
  --ml Minimum length   Minimum length of reads kept after merging [default:
                        380]

```

## Examples:


##### Using silva database:
```
$ auto-q.py -i /data/experiment1/fastqs/ -o /data/experiment1/results/ -t 12 -p 10 -r silva -n 10 -e 5000 -c /bin/auto-q/qiime.cfg 
```

##### Stop at merging step:
```buildoutcfg
$ auto-q.py -i /data/experiment1/fastqs/ -o /data/experiment1/results/ -t 10 -p 16 -s merging -n 10 -c /bin/auto-q/qiime.cfg
```

##### Begin analysis with otu picking using fasta files after chimera removal:
```buildoutcfg
$ auto-q.py -i /data/experiment1/results/chi/ -o /data/experiment1/results/ -b otu_picking -n 10 -c /bin/auto-q/qiime.cfg
```

##### If not going run the script in parallel use (-n 1):

```buildoutcfg
$ auto-q.py -i /data/experiment1/results/chi/ -o /data/experiment1/results/ -b otu_picking -n 1 -c /bin/auto-q/qiime.cfg
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

## Stop at:



## How to cite:
The paper of the script is under preparation now.
