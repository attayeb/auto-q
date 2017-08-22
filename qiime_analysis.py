#! /usr/bin/env python
from __future__ import unicode_literals
# coding: utf-8
# import sys

## more to do:
# 1. log file
# 2. save parameter file and configuration file in results folder
#

## updates:
# 0.2.4 Using fastq-join with trimming 10 and p 16 with silva as default
import logging
import datetime
import os  # operation system packages
import configparser  # a package to parse INI file or confige file.
import argparse  # a package to parse commandline arguments.
import pandas as pd  # pandas package

from subprocess import call  # to run command line scripts
from subprocess import Popen, PIPE
#from pathos.multiprocessing import ProcessingPool as Pool  # for parallel computing
from multiprocessing.dummy import Pool as Pool

__version__ = '0.2.4'
__author__ = "Attayeb Mohsen"
__date__ = "18/8/2017"

starting_message = """ Microbiome analysis using multiple methods
    Version: %s
    Date: %s 
    National Institutes of Biomedical Innovation, Helath, and Nutrition\n"""%(__version__, __date__)

ID = str(datetime.datetime.now())
ID = ID.replace(" ", "")
ID = ID.replace(":", "")
ID = ID.replace(".", "")
ID = ID.replace("-", "")
ID = ID[0:14]

PR = dict({"id": ID})  # PARAMETERS dict


def asfolder(folder):
    if folder[-1] != "/":
        return (folder + "/")
    else:
        return (folder)


def execute(command, shell):
    p = Popen(command.split(), stderr=PIPE, stdout=PIPE)
    output, error = p.communicate()
    loginfo(command)
    if output != b"":
        loginfo(output.encode('utf-8'))
    if error != b"":
        logwarning(error.encode('utf-8'))


def loginfo(message):
    logging.info(message.encode('utf-8'))


def logwarning(message):
    logging.warning(message.encode('utf-8'))


def get_configuration():
    global PR
    cp = configparser.ConfigParser()
    cp.read(PR['ConfigFile'])
    PR['Ftrimmed'] = asfolder(cp.get('FOLDERS', 'trimmed'))
    PR['Fmerged'] = asfolder(cp.get('FOLDERS', 'merged'))
    PR['Fqc'] = asfolder(cp.get('FOLDERS', 'quality_step'))
    PR['Fchi'] = asfolder(cp.get('FOLDERS', 'chimera_removed'))
    PR['Fotus'] = asfolder(cp.get('FOLDERS', 'otus'))
    PR['Fdiv'] = asfolder(cp.get('FOLDERS', 'diversity_analyses'))
    PR['Fothers'] = asfolder(cp.get('FOLDERS', 'others'))
    PR['number_of_cores'] = int(cp.get('GENERAL', 'jobs_to_start'))
    PR['silva_taxonomy'] = cp.get('SILVA', 'taxonomy')
    PR['silva_reference_seqs'] = cp.get('SILVA', 'reference_seqs')
    PR['silva_core_alignment'] = cp.get('SILVA', 'core_alignment')
    PR['silva_chim_ref'] = cp.get('CHIMERA', 'silva')
    PR['gg_chim_ref'] = cp.get('CHIMERA', 'gg')
    PR['unite_taxonomy'] = cp.get('UNITE', 'taxonomy')
    PR['unite_reference_seqs'] = cp.get('UNITE', 'reference_seqs')
    PR['similarity'] = cp.get('GENERAL', 'similarity')
    PR['blast_e_value'] = cp.get('GENERAL', 'blast_e_value')



def write_parameter_file(parameter_file):
    if PR['silva'] == True:
        parameter_string = """
    assign_taxonomy:id_to_taxonomy_fp\t%(taxonomy)s
    assign_taxonomy:reference_seqs_fp\t%(reference_seqs)s
    pick_otus.py:pick_otus_reference_seqs_fp\t%(reference_seqs)s
    filter_alignment.py:pynast_template_alignment_fp\t%(core_alignment)s
    parallel:jobs_to_start\t%(jobs_to_start)d
    assign_taxonomy:similarity\t%(similarity)s
    """ % {'taxonomy': PR['silva_taxonomy'],
           'reference_seqs': PR['silva_reference_seqs'],
           'core_alignment': PR['silva_core_alignment'],
           'jobs_to_start': PR['number_of_cores'],
           'similarity': PR['similarity']}
    elif PR['fungus'] == True:
        # pass
        parameter_string = """
        assign_taxonomy:id_to_taxonomy_fp\t%(taxonomy)s
        assign_taxonomy:reference_seqs_fp\t%(reference_seqs)s
        pick_otus.py:pick_otus_reference_seqs_fp\t%(reference_seqs)s
        parallel:jobs_to_start\t%(jobs_to_start)d
        assign_taxonomy:assignment_method blast
        # should we use e_value or blast_e_value
        parallel_assign_taxonomy_blast:e_value\t%(blast_e_value)s
        # comment
        """ % {'taxonomy': PR['unite_taxonomy'],
               'reference_seqs': PR['unite_reference_seqs'],
               'jobs_to_start': PR['number_of_cores'],
               'blast_e_value': PR['blast_e_value']}
    else:
        parameter_string = '''
    parallel:jobs_to_start\t%(jobs_to_start)d
            ''' % {'jobs_to_start': PR['number_of_cores']}
    call("mkdir -p %s" % PR['others'] , shell=True)
    f = open(parameter_file, "w")
    f.write(parameter_string)
    f.close()


def trimfolder(in_folder, out_folder, trimq):
    """
    Trim the right side of fastq files using BBDUK,
    in_folder: input is a folder of fastq files with pair ends
    both pairs will be trimmed together
    out_folder is the folder were trimmed files is to be saved
    trimq is the threshold of trimming: 12 is recommended
    trimfolder(in_folder, out_folder, trimq)
    ex:
    trimfolder("/experiment1/data/fastq/",
           "/experiment1/analysis/trimmed/", trimq=12)

    """

    import os
    ## step to ensure that "/" is at the end of the folders'
    ## names
    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    files = os.listdir(in_folder)
    files.sort()
    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]
    call("mkdir -p %s" % out_folder, shell=True)

    # get_ipython().system(u'mkdir -p {out_folder}')
    def process(i):

        in1 = in_folder + ins1[i]
        in2 = in_folder + ins2[i]
        print("Trimming: %s and %s" % (ins1[i], ins2[i]))
        out1 = out_folder + ins1[i]
        out2 = out_folder + ins2[i]
        execute(
            "bbduk.sh -Xmx1000m -in1=%s -in2=%s -out1=%s -out2=%s -qtrim=r -trimq=%d" % (in1, in2, out1, out2, trimq),
            shell=True)

    p = Pool(number_of_cores)
    p.map(process, range(len(ins1)))


def mergefolderbb(in_folder, out_folder, maxloose=True):
    """
    Merge a folder of pair ends fastq files using BBMERGE
    in_folder : input
    out_folder : output
    maxloose is set as default, if you want to use
    default parameters assign maxloose=False

    mergefolderbb(in_folder, out_folder, maxloose=True)

    ex:
    maxloose:
    mergefolderbb("/data/fastq/", "data/merged/", maxloose=True)
    default:
    mergefolderbb("/data/fastq/", "data/merged/", maxloose=False)


    """
    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    import os
    files = os.listdir(in_folder)
    files.sort()

    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]
    outs = [x.replace("_L001_R1_001", "") for x in ins1]
    # get_ipython().system(u'mkdir -p {out_folder}')
    call("mkdir -p %s" % out_folder, shell=True)

    def process(i):

        in1 = in_folder + ins1[i]
        in2 = in_folder + ins2[i]
        print("Merging: %s and %s" % (ins1[i], ins2[i]))
        out = out_folder + outs[i]
        if maxloose:
            # get_ipython().system(
            #    u'bbmerge.sh -in1={in1} -in2={in2} -out={out} -maxloose=t -ignorebadquality')
            execute("bbmerge.sh -in1=%s -in2=%s -out=%s -maxloose=t -ignorebadquality" % (in1, in2, out), shell=True)



        else:
            # get_ipython().system(u'bbmerge.sh -in1={in1} -in2={in2} -out={out} -ignorebadquality')
            execute("bbmerge.sh -in1=%s -in2=%s -out=%s -ignorebadquality" % (in1, in2, out), shell=True)

    p = Pool(number_of_cores)
    p.map(process, range(len(ins1)))
    print("Merging finished.")

def mergefolderfastq(in_folder, out_folder, pp):
    """
    Merge folder of fastq files following illumina naming system
    using fastq-join
    mergefolderfastq(in_folder, out_folder, pp)
    in_folder: the source folder of fastq files.
    out_folder: the target folder for merged files.
    p: percentage of maximum difference.
    ex:
        mergefolderfastq("/home/data/human/",
        "/home/data/merged-human/", p=12)
    """

    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    files = os.listdir(in_folder)

    files.sort()

    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]

    # modify the output file names
    outs = [x.replace("_L001_R1_001", "") for x in ins1]
    # get_ipython().system(u'mkdir -p {out_folder}')
    call("mkdir -p %s" % out_folder, shell=True)

    # for i in range(len(ins1)):
    def process(i):
        in1 = in_folder + ins1[i]
        in2 = in_folder + ins2[i]
        print("Merging: %s and %s " % (ins1[i], ins2[i]))
        out = out_folder + outs[i]
        # get_ipython().system(u'fastq-join -p {pp} {in1} {in2} -o {out}')
        execute("fastq-join -p %d %s %s -o %s" % (pp, in1, in2, out), shell=True)

        # get_ipython().system(u'rm {out+"un1"}')
        call("rm %sun1" % out, shell=True)
        # get_ipython().system(u'rm {out+"un2"}')
        call("rm %sun2" % out, shell=True)
        # get_ipython().system(u'mv {out+"join"} {out}')
        call("mv %sjoin %s" % (out, out), shell=True)

    print("Merging finished")

    p = Pool(number_of_cores)
    p.map(process, range(len(ins1)))


def qualitycontrol(in_folder, out_folder, q):
    """
    Using split_libraries script from QIIME
    Out_put is fasta files
    all other files will be suppressed
    q: is the quality threshold

    qualitycontrol(in_folder, out_folder, q)
    qualitycontrol("data1/merged/", "/data1/qc/", q=19)

    """
    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    import os
    files = os.listdir(in_folder)
    files.sort()

    # get_ipython().system(u'mkdir -p {temp}')
    call("mkdir -p %s " % out_folder, shell=True)

    def process(i):
        # for i in files:
        temp = out_folder + "temp" + i + "/"
        print("Quality control: %s" % i)
        sampleid = i.replace(".fastq", "")
        infile = in_folder + i
        outfile = out_folder + i.replace(".fastq", ".fasta")
        # get_ipython().system(
        #    u"split_libraries_fastq.py -i {infile} -o {temp}         --barcode_type 'not-barcoded' --sample_ids {sampleid} -q {q}")
        execute("""split_libraries_fastq.py -i %s -o %s --barcode_type not-barcoded --sample_ids %s -q %s""" % (
            infile, temp, sampleid, q), shell=True)

        tempfile = temp + "seqs.fna"
        # get_ipython().system(u'mv {tempfile} {outfile}')
        call("mv %s %s" % (tempfile, outfile), shell=True)
        call("rm -r %s" % temp, shell=True)

    p = Pool(5)
    p.map(process, files)
    print("Quality control finished.")

    # get_ipython().system(u'rm -r {temp}')


def removechimera(in_folder, out_folder, silva=False):
    """
    removechimera(in_folder,out_folder, silva=False)
    GreenGenes database are used as default
    if you want to use silve assign silva=True

    ex:
    using Silva:
       removechimera("/data/qc/", "/data/chi/", silva=True)
    using GG:
       removechimera("/data/qc/", "/data/chi/", silva=False)
       or:
       removechimera("/data/qc/", "/data/chi/")

    """
    global PR
    import os
    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    files = os.listdir(in_folder)
    files.sort()

    # get_ipython().system(u'mkdir -p {temp}')
    call("mkdir -p %s" % out_folder, shell=True)

    def process(i):
        print("Chimera removal: %s" %i)
        temp = out_folder + "temp" + i + "/"
        if silva:
            execute("identify_chimeric_seqs.py -i %s -m usearch61 -o %s -r %s"
                    % (in_folder + i, temp + i, PR['silva_chim_ref']),
                    shell=True)
        else:
            execute("identify_chimeric_seqs.py -i %s -m usearch61 -o %s -r %s" % (
                in_folder + i, temp + i, PR['gg_chim_ref']),
                    shell=True)

            # get_ipython().system(
        # u'filter_fasta.py -f {in_folder}{i} -o {out_folder}{i}         -s {temp}{i}/non_chimeras.txt')
        execute("filter_fasta.py -f %s -o %s -s %s/non_chimeras.txt" % (in_folder + i, out_folder + i, temp + i),
                shell=True)
        call("rm -r %s" % temp, shell=True)

    # get_ipython().system(u'rm -r {temp}')
    p = Pool(number_of_cores)
    p.map(process, files)


def pickotus(in_folder, out_folder, parameter_file, silva=False, fungus=False):
    """
    pickotus(in_folder, out_folder, silva=False)
    input (in_folder): is a folder of fasta files
    if you want to use silve assign it as True
    ex:
    pickotus(in_folder, out_folder, silva=True)

    """

    global PR
    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    infolder = in_folder + "*.fasta"
    print("Otu picking...")
    if silva:

        # get_ipython().system(
        #    u'pick_open_reference_otus.py -i "$infolder"         -o {out_folder} -p otu_pick_param         -r /home/qiime/Database/SILVA_128_QIIME_release/rep_set/rep_set_16S_only/97/97_otus_16S.fasta')
        execute("pick_open_reference_otus.py -i %s -o %s -p %s -r %s -a -O %d "
                % (infolder, out_folder, PR['parameter_file_name'], PR['silva_reference_seqs'], PR['number_of_cores']),
                shell=True)


    elif fungus:
        execute("pick_open_reference_otus.py -i %s -o %s -p %s -a -O %d --supress_align_and_tree"
                % (infolder, out_folder, PR['parameter_file_name'], PR['number_of_cores']), shell=True)

    else:
        # get_ipython().system(u'pick_open_reference_otus.py -i "$infolder" -o {out_folder} -p otu_para_gg.txt')
        execute("pick_open_reference_otus.py -i %s -o %s -p %s -a -O %d"
                % (infolder, out_folder, PR['parameter_file_name'], PR['number_of_cores']), shell=True)


def create_map(in_folder, out_file):
    """
    create_map(in_folder, out_file)
    in_folder: is the folder: results from otu picking
    out_file is the file name !! Not folder
    Create basic mapping file,
    If you have other parameters, please use excell or
    libreoffice to add them, like for example Age, or Treated

    ex:
    create_map("/human/experiment22/analyzed/otus/",
    "/human/experiments/analyzed/samples_map.txt")

    """
    in_folder = asfolder(in_folder)
    print("Writing mapping file")
    import os
    sampleids = os.listdir(in_folder)
    ids = [x.replace(".fasta", "") for x in sampleids]
    ids = [x.split("_")[0] for x in ids]
    d = {'#SampleID': ids}
    df = pd.DataFrame(data=d)
    df['BarcodeSequence'] = ""
    df['LinkerPrimerSequence'] = ""
    df['Read'] = "R1"
    df["File"] = sampleids
    df["Description"] = "single_file"
    df.to_csv(out_file, sep=str(u"\t"), index=False)


def corediv(in_folder, out_folder, mapping_file, depth):
    """
    corediv(in_folder, out_folder, mapping_file, depth)
    in_folder: the otus folder
    depth of subsampling

    ex:
    corediv("/human/aging/otus/", "/human/aging/div/",
    "/human/aging/sample_map.txt", 10000)
    """
    print("Core diversity analyses...")
    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)
    if PR['fungus']:
        biom = in_folder + "otu_table_mc2_w_tax.biom"
    else:
        biom = in_folder + "otu_table_mc2_w_tax_no_pynast_failures.biom"

    tree = in_folder + "rep_set.tre"
    # get_ipython().system(
    #    u'core_diversity_analyses.py -i {biom}     -o {out_folder}     -m {mapping_file}     -t {tree}     -e {depth}')
    if PR['fungus']:
        execute("core_diversity_analyses.py -i %s -o %s -m %s -e %d --nonphylogenetic_diversity" % (biom, out_folder, mapping_file, depth),
            shell=True)
    else:
        execute(
            "core_diversity_analyses.py -i %s -o %s -m %s -t %s -e %d" % (biom, out_folder, mapping_file, tree, depth),
            shell=True)


def full_analysis(in_folder, out_folder, depth, silva, trimq, fastq_join,
                  qcq, maxloose, fastq_p):
    global PR
    """
    Full analysis will be done
    using bbmerge with maxloose
    quality control with q=19
    using gg database for otu-picking and chimera removal
    set the depth of diversity analysis to 10000
    """

    trimmed = asfolder(out_folder + PR['Ftrimmed'])
    merged = asfolder(out_folder + PR['Fmerged'])
    qc = asfolder(out_folder + PR['Fqc'])
    chi = asfolder(out_folder + PR['Fchi'])
    otus = asfolder(out_folder + PR['Fotus'])
    div = asfolder(out_folder + PR['Fdiv'])

    trimfolder(in_folder, trimmed, trimq)
    if fastq_join:
        mergefolderfastq(trimmed, merged, fastq_p)
    else:
        mergefolderbb(trimmed, merged, maxloose=maxloose)

    qualitycontrol(merged, qc, qcq)
    removechimera(qc, chi, silva)
    pickotus(chi, otus, silva)
    create_map(qc, PR['others'] + "map.tsv")
    corediv(otus, div, PR['others'] + "map.tsv", depth)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Microbiome analysis using multiple methods
    Version: %s
    Date: %s """%(__version__, __date__))
    parser.add_argument("-i",  # "--input",
                        dest="input",
                        # type=str,
                        help="The folder where Fastq files are stored [required]",
                        required=True)

    parser.add_argument("-o",
                        # "--output",
                        dest="output",
                        type=str,
                        help="The folder of all results [required]",
                        required=True)

    # parser.add_argument("-j",
    #                    # "--output",
    #                    dest="job_name",
    #                    type=str,
    #                    help="The name of the job")


    parser.add_argument("-t",
                        # "--trim_phred_quality_threshold",
                        dest="trim_threshold",
                        type=int,
                        help="phred quality threshold for trimming [10]",
                        default=10)

    parser.add_argument("-j",
                        # "--fastq-join",
                        dest='fastq_join',
                        help="use fastq-join for joining [False]",
                        action="store_true",
                        default=True)

    parser.add_argument("-d",
                        # "--percentage_of_maximum_difference",
                        type=int,
                        dest="fastq_p",
                        help="Percentage of maximum difference in fastq-join [16]",
                        default=16)


    parser.add_argument("-q",
                        # "--quality_control_phred_threshold",
                        dest="qc_threshold",
                        type=int,
                        help="quality control phred threshold [19]",
                        default=19)

    # parameter file

    parser.add_argument("-c",
                        # "--ConfigFile",
                        dest="ConfigFile",
                        type=str,
                        default='qiime.cfg',
                        help="Configuration file name [qiime.cfg]")

    # parser.add_argument("-a",
    #                    "--automatic_parameter_file",
    #                    dest = "automatic_parameter_file",
    #                    help="Automatically produce parameter file \
    #                    using information from configuration file",
    #                    action="store_true")

    parser.add_argument("-p",
                        # "--parameter_file_name",
                        help="The name of the parameter file [if not assigned is automatically produced using "
                             "configuration file",
                        type=str,
                        dest="parameter_file_name")

    parser.add_argument("-n",
                        # "--number_of_cores",
                        help="Number of cores to be used for the analysis [2]",
                        type=int,
                        dest="number_of_cores",
                        default=2)


    parser.add_argument("-f",
                        # "--silva",
                        dest="fungus",
                        help="Using Unite database for fungal samples[False]",
                        action="store_true")

    parser.add_argument("-m",
                        # "--maxloose",
                        dest="maxloose",
                        help="Assign maxloose to be true for bbmerge [False]",
                        action="store_true")


    parser.add_argument("-s",
                        # "--silva",
                        dest="silva",
                        help="Using SILVA database [False]",
                        action="store_true",
                        default=True)


    parser.add_argument("-e",
                        # "--depth",
                        dest="depth",
                        type=int,
                        help="set the depth of diversity analyses [10000]",
                        default=10000)



    x = parser.format_usage()
    parser.usage = starting_message + x
    arg = parser.parse_args()

    # Namespace(fastq_join=False, fastq_p=None, input='/h', maxloose=False, output='/e', qc_threshold=None, silva=False, trim_threshold=None)
    ## global parameters assignment::


    PR.update({
        'in_folder': arg.input,
        'out_folder': arg.output,
        'silva': arg.silva,
        'qcq': arg.qc_threshold,
        'maxloose': arg.maxloose,
        'trimq': arg.trim_threshold,
        'fastq_join': arg.fastq_join,
        'fastq_p': arg.fastq_p,
        'depth': arg.depth,
        'fungus': arg.fungus,
        'ConfigFile': arg.ConfigFile,
        # 'automatic_parameter_file': arg.automatic_parameter_file,
        'parameter_file_name': arg.parameter_file_name})

    ## parameter_file
    get_configuration()
    PR['others'] = asfolder(PR['out_folder'] + PR['Fothers'])
    if arg.parameter_file_name == None:
        PR['parameter_file_name'] = PR['others']+"para%s.txt" % PR['id']
        write_parameter_file(PR['parameter_file_name'])
    logging.basicConfig(filename=PR['others']+"log.txt",
                        format='%(asctime)s %(levelname)s \n %(message)s',
                        level = logging.DEBUG)

    number_of_cores = PR['number_of_cores']
    loginfo('started')
    [loginfo(str(P)+": "+str(PR[P])) for P in PR]



    full_analysis(in_folder=PR['in_folder'],
                  out_folder=PR['out_folder'],
                  silva=PR['silva'],
                  fastq_join=PR['fastq_join'],
                  fastq_p=PR['fastq_p'],
                  maxloose=PR['maxloose'],
                  qcq=PR['qcq'],
                  depth=PR['depth'],
                  trimq=PR['trimq'])

    loginfo("Finished")
