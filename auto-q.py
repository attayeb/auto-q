#! /usr/bin/env python
from __future__ import unicode_literals
# coding: utf-8
# import sys

## more to do:


# 1. log file
# 2. save parameter file and configuration file in results folder
# 3. merge fasta files before closed reference otu picking or use filtering
## 4. remove primers from 3' end and 5'end (length will be decided by Park san)
## 5. remove short reads < 380
## remove phix should be done using bbduk with phix174_ill.ref.fadrr


## updates:
# 0.2.4 Using fastq-join with trimming 10 and p 16 with silva as default
## 0.2.5: added close reference otu picking using -reference tag: default is open # removed
## 0.2.5: added beginwith: start with either otu_picking or diversity_analysis
## with diversity analysis mapping file is required
## 0.2.6 added force trim left for R1 .. 21 bases in trimming step using python script
## 0.2.6 added force trim left for R2 .. 17 bases in trimming step using python script
## 0.2.6 remove short reads < 380
## 0.2.6 remove adapter sequence using bbduk



import logging
import datetime
import os               # operation system packages
import configparser     # a package to parse INI file or confige file.
import argparse         # a package to parse commandline arguments.
import sys
from re import sub

from subprocess import call  # to run command line scripts
from subprocess import Popen, PIPE, check_output
from multiprocessing.dummy import Pool as Pool

__version__ = '0.2.6'
__author__ = "Attayeb Mohsen"
__date__ = "19/10/2017"

starting_message = """ Microbiome analysis using multiple methods
    Version: %s
    Date: %s 
    National Institutes of Biomedical Innovation, Health, and Nutrition\n""" % (__version__, __date__)

ID = str(datetime.datetime.now())
ID = ID.replace(" ", "")
ID = ID.replace(":", "")
ID = ID.replace(".", "")
ID = ID.replace("-", "")
ID = ID[0:14]

PR = dict({"id": ID})  # PARAMETERS dict


def remove_short_reads(infqfile, outfqfile, length):
    """

    :param infqfile: input fastq file name.
    :type infqfile: str
    :param outfqfile: output fastq file name, after removing short reads.
    :type outfqfile: str
    :param length: minimum reads length.
    :type length: int
    :rtype: None
    :return: None
    @Action: filter fastq files removing short reads

    """
    infq = open(infqfile, "r")
    outfq = open(outfqfile, "w")
    lines = infq.readlines()
    for a, b, c, d in zip(lines[0::4], lines[1::4], lines[2::4], lines[3::4]):
        if len(b) > length:
            outfq.write(a)
            outfq.write(b)
            outfq.write(c)
            outfq.write(d)

    infq.close()
    outfq.close()

def asfolder(folder):
    """
    Add "/" at the end of the folder if not inserted

    :param folder: the folder name
    :type folder: str
    :return: file names with / at the end
    :rtype: str
    """
    if folder[-1] != "/":
        return (folder + "/")
    else:
        return (folder)



def execute(command, shell):
    """
    Execute command using os package and return output to log file

    :param command: The command to be executed
    :type command: str
    :param shell: Takes either True or False
    :type shell: boolean
    :return: Run the command in the background and save the output to the logging file.

    """

    p = Popen(command.split(), stderr=PIPE, stdout=PIPE)
    output, error = p.communicate()
    loginfo(command)
    if output != b"":
        loginfo(output.encode('utf-8'))
    if error != b"":
        logwarning(error.encode('utf-8'))



def loginfo(message):
    """
    save information to log file

    :param message: saved to log file
    :type message: str
    :return:
    """
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
    PR['gg_taxonomy'] = cp.get('GG', 'taxonomy')
    PR['gg_reference_seqs'] = cp.get('GG', 'reference_seqs')
    PR['gg_core_alignment'] = cp.get('GG', 'core_alignment')
    PR['gg_chim_ref'] = cp.get('CHIMERA', 'gg')
    PR['unite_taxonomy'] = cp.get('UNITE', 'taxonomy')
    PR['unite_reference_seqs'] = cp.get('UNITE', 'reference_seqs')
    PR['similarity'] = cp.get('GENERAL', 'similarity')
    PR['blast_e_value'] = cp.get('GENERAL', 'blast_e_value')
    PR['bbmap_resources'] = cp.get('bbmap', 'resources')


def locate_bbmap():
    """

    :return:
    """
    folder = check_output(["locate", "bbmerge.sh"]).decode("utf-8")
    return (sub('bbmerge.sh\n$', '', folder))


def check_before_start():
    """

    :return:
    """
    if os.path.isfile(PR['ConfigFile']):
        pass
    else:
        raise IOError("configuration file does not exist")

    if PR['rdb'] == "silva":
        condition = True
        condition = condition and os.path.isfile(PR['silva_taxonomy'])
        condition = condition and os.path.isfile(PR['silva_reference_seqs'])
        condition = condition and os.path.isfile(PR['silva_core_alignment'])
        if not condition:
            raise IOError("Can not find Silva database files, "
                          "please check the configuration file: %s to set up the correct folder" % PR['ConfigFile'])
    if PR['rdb'] == "gg":
        condition = True
        condition = condition and os.path.isfile(PR['gg_taxonomy'])
        condition = condition and os.path.isfile(PR['gg_reference_seqs'])
        condition = condition and os.path.isfile(PR['gg_core_alignment'])
        if not condition:
            raise IOError("Can not find greengenes database files, "
                          "please check the configuration file: %s to set up the correct folder" % PR['ConfigFile'])
    if os.path.isdir(PR['out_folder']):
        raise IOError("Output folder exists, Please use a non existent folder name")


def write_parameter_file(parameter_file):
    """

    :param parameter_file:
    :return:
    """
    if PR['rdb'] == "silva":
        parameter_string = """
    assign_taxonomy:id_to_taxonomy_fp\t%(taxonomy)s
    assign_taxonomy:reference_seqs_fp\t%(reference_seqs)s
    pick_otus.py:pick_otus_reference_seqs_fp\t%(reference_seqs)s
    pick_otus:enable_rev_strand_match True
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
    assign_taxonomy:id_to_taxonomy_fp\t%(taxonomy)s
    assign_taxonomy:reference_seqs_fp\t%(reference_seqs)s
    pick_otus.py:pick_otus_reference_seqs_fp\t%(reference_seqs)s
    pick_otus:enable_rev_strand_match True
    filter_alignment.py:pynast_template_alignment_fp\t%(core_alignment)s
    parallel:jobs_to_start\t%(jobs_to_start)d
    assign_taxonomy:similarity\t%(similarity)s
    ''' % {'taxonomy': PR['gg_taxonomy'],
           'reference_seqs': PR['gg_reference_seqs'],
           'core_alignment': PR['gg_core_alignment'],
           'jobs_to_start': PR['number_of_cores'],
           'similarity': PR['similarity']}
    os.mkdir(PR['others'])

    f = open(parameter_file, "w")
    f.write(parameter_string)
    f.close()


def primertrim(infqfile, outfqfile, length):
    """

    :param infqfile:
    :param outfqfile:
    :param length:
    :return:
    """

    infq = open(infqfile, "r")
    outfq = open(outfqfile, "w")
    lines = infq.readlines()
    for a, b, c, d in zip(lines[0::4], lines[1::4], lines[2::4], lines[3::4]):
        outfq.write(a)
        outfq.write(b[length:])
        outfq.write(c)
        outfq.write(d[length:])

    infq.close()
    outfq.close()


def trimfolder(in_folder, out_folder, trimq, ftrim=True):
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

    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    files = os.listdir(in_folder)
    files.sort()
    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]
    os.mkdir(out_folder)
    # call("mkdir -p %s" % out_folder, shell=True)
    print("Trimming...")

    # get_ipython().system(u'mkdir -p {out_folder}')
    def process(i):
        in1 = in_folder + ins1[i]
        in2 = in_folder + ins2[i]
        print("%s and %s" % (ins1[i], ins2[i]))
        out1 = out_folder + ins1[i]
        out2 = out_folder + ins2[i]
        out1_temp1 = out_folder + "temp1_" + ins1[i]
        out2_temp1 = out_folder + "temp1_" + ins2[i]

        # forctrimleft was added
        if ftrim:
            primertrim(in1, out1_temp1, 17)
            primertrim(in2, out2_temp1, 21)
            # execute(
            #    "bbduk.sh -Xmx1000m -in=%s -out=%s -forcetrimleft=17" % (in1, out1_temp1),
            #    shell=True)
            # execute(
            #    "bbduk.sh -Xmx1000m -in=%s -out=%s -forcetrimleft=21" % (in2, out2_temp1),
            #    shell=True)
            #execute(
            #    "bbduk.sh -Xmx1000m -in1=%s -in2=%s -out1=%s -out2=%s -qtrim=r -trimq=%d" %
            #    (out1_temp1, out2_temp1, out1_temp2, out2_temp2, trimq),
            #    shell=True)
            #os.remove(out1_temp1)
            #os.remove(out2_temp1)
        else:
            out1_temp1 = in1
            out2_temp1 = in2

        if PR['adapter_ref'] != None:
            #print(PR['adapter_ref'])
            #execute(
            #    "bbduk.sh -Xmx1000m -in1=%s -in2=%s -out1=%s -out2=%s ref=%s, ktrim=r k=23 mink=11 hdist=1 tpe tbo" %
            #    (out1_temp2, out2_temp2, out1, out2, PR['adapter_ref']), shell=True)

            execute("bbduk.sh -Xmx1000m -in1=%s -in2=%s -out1=%s -out2=%s -outm=stdout.fa -ref=%s -qtrim=r -trimq=%d -k=18 -ktrim=f" %
                (out1_temp1, out2_temp1, out1, out2, PR['adapter_ref'], trimq), shell=True)
        else:
            execute(
                "bbduk.sh -Xmx1000m -in1=%s -in2=%s -out1=%s -out2=%s -qtrim=r -trimq=%d" %
                (out1_temp1, out2_temp1, out1, out2, trimq), shell=True)

        os.remove(out1_temp1)
        os.remove(out2_temp1)

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
    os.mkdir(out_folder)
    # call("mkdir -p %s" % out_folder, shell=True)
    print("Merging ...")

    def process(i):

        in1 = in_folder + ins1[i]
        in2 = in_folder + ins2[i]


        print("%s and %s" % (ins1[i], ins2[i]))
        out = out_folder + outs[i]
        if maxloose:
            execute("bbmerge.sh -in1=%s -in2=%s -out=%s -maxloose=t -ignorebadquality" % (in1, in2, out), shell=True)



        else:
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
    global PR
    in_folder = asfolder(in_folder)
    out_folder = asfolder(out_folder)

    files = os.listdir(in_folder)

    files.sort()

    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]

    outs = [x.replace("_L001_R1_001", "") for x in ins1]
    os.mkdir(out_folder)

    # call("mkdir -p %s" % out_folder, shell=True)


    def process(i):
        in1 = in_folder + ins1[i]
        in2 = in_folder + ins2[i]
        print("Merging: %s and %s " % (ins1[i], ins2[i]))
        out = out_folder + "temp_" + outs[i]
        out_final = out_folder + outs[i]
        # out = out_final
        execute("fastq-join -p %d %s %s -o %s" % (pp, in1, in2, out), shell=True)
        os.remove("%sun1" % out)
        os.remove("%sun2" % out)
        os.rename("%sjoin" % out, out)
        remove_short_reads(out, out_final, PR['minimum_length'])
        os.remove(out)

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
    os.mkdir(out_folder)

    # call("mkdir -p %s " % out_folder, shell=True)

    def process(i):
        temp = out_folder + "temp" + i + "/"
        print("Quality control: %s" % i)
        sampleid = i.replace(".fastq", "")
        infile = in_folder + i
        outfile = out_folder + i.replace(".fastq", ".fasta")
        execute("""split_libraries_fastq.py -i %s -o %s --barcode_type not-barcoded --sample_ids %s -q %s""" % (
            infile, temp, sampleid, q), shell=True)

        tempfile = temp + "seqs.fna"
        call("mv %s %s" % (tempfile, outfile), shell=True)
        call("rm -r %s" % temp, shell=True)

    p = Pool(5)
    p.map(process, files)
    print("Quality control finished.")


def removechimera(in_folder, out_folder, rdb="silva"):
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

    os.mkdir(out_folder)

    # call("mkdir -p %s" % out_folder, shell=True)

    def process(i):
        print("Chimera removal: %s" % i)
        temp = out_folder + "temp" + i + "/"
        if rdb == "silva":
            execute("identify_chimeric_seqs.py -i %s -m usearch61 -o %s -r %s"
                    % (in_folder + i, temp + i, PR['silva_chim_ref']),
                    shell=True)
        else:
            execute("identify_chimeric_seqs.py -i %s -m usearch61 -o %s -r %s" % (
                in_folder + i, temp + i, PR['gg_chim_ref']),
                    shell=True)

        execute("filter_fasta.py -f %s -o %s -s %s/non_chimeras.txt" % (in_folder + i, out_folder + i, temp + i),
                shell=True)
        call("rm -r %s" % temp, shell=True)

    p = Pool(number_of_cores)
    p.map(process, files)


def pickotus(in_folder, out_folder, rdb="silva", fungus=False):
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

    if rdb == "silva":
        execute("pick_open_reference_otus.py -i %s -o %s -p %s -r %s -a -O %d "
                % (infolder, out_folder, PR['parameter_file_name'], PR['silva_reference_seqs'], PR['number_of_cores']),
                shell=True)
        execute("filter_otus_from_otu_table.py -i %s -o %s --negate_ids_to_exclude -e %s"
                % (out_folder + "otu_table_mc2_w_tax_no_pynast_failures.biom",
                   out_folder + "otu_table_mc2_w_tax_no_pynast_failures_close_reference.biom",
                   PR['silva_reference_seqs']), shell=True)

    elif fungus:
        execute("pick_open_reference_otus.py -i %s -o %s -p %s -a -O %d --supress_align_and_tree"
                % (infolder, out_folder, PR['parameter_file_name'], PR['number_of_cores']), shell=True)

    else:
        execute("pick_open_reference_otus.py -i %s -o %s -r %s -p %s -a -O %d"
                % (infolder, out_folder,
                   PR['gg_reference_seqs'], PR['parameter_file_name'],
                   PR['number_of_cores']), shell=True)

        execute("filter_otus_from_otu_table.py -i %s -o %s --negate_ids_to_exclude -e %s"
                % (out_folder + "otu_table_mc2_w_tax_no_pynast_failures.biom",
                   out_folder + "otu_table_mc2_w_tax_no_pynast_failures_close_reference.biom",
                   PR['gg_reference_seqs']), shell=True)


def writedf(out_file, ids, sampleids):
    f = open(out_file, "w+")
    f.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tRead\tFile\tDescription\n")
    for x in range(len(ids)):
        f.write("%s\t\t\tR1\t%s\tsingle_file\n"%(ids[x],sampleids[x] ))
    f.close()




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
    writedf(out_file, ids, sampleids)


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
        biom_singleton = in_folder + "otu_table_mc2_w_tax_no_pynast_failures_singleton_filtered.biom"
    tree = in_folder + "rep_set.tre"
    # get_ipython().system(
    #    u'core_diversity_analyses.py -i {biom}     -o {out_folder}     -m {mapping_file}     -t {tree}     -e {depth}')
    if PR['fungus']:
        execute("core_diversity_analyses.py -i %s -o %s -m %s -e %d --nonphylogenetic_diversity" % (
            biom, out_folder, mapping_file, depth),
                shell=True)
    else:
        execute(
            "core_diversity_analyses.py -i %s -o %s -m %s -t %s -e %d" % (biom, out_folder, mapping_file, tree, depth),
            shell=True)


def full_analysis(in_folder, out_folder, depth, rdb, trimq, joining_method,
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
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("error method")
    qualitycontrol(merged, qc, qcq)
    removechimera(qc, chi, rdb)
    pickotus(chi, otus, rdb)
    # here
    if create_mapping_file:
        create_map(qc, PR['mapping_file'])
    corediv(otus, div, PR['mapping_file'], depth)


def stop_at_merging(in_folder, out_folder, trimq, joining_method, maxloose, fastq_p):
    global PR
    trimmed = asfolder(out_folder + PR['Ftrimmed'])
    merged = asfolder(out_folder) + PR['Fmerged']
    trimfolder(in_folder, trimmed, trimq)
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("%s: unknown merging metod method" % joining_method)


def stop_at_quality_control(in_folder, out_folder, joining_method, trimq,
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

    trimfolder(in_folder, trimmed, trimq)
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("%s: unknown merging metod method" % joining_method)
    qualitycontrol(merged, qc, qcq)


def stop_at_chimera_removal(in_folder, out_folder, rdb, trimq, joining_method,
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

    trimfolder(in_folder, trimmed, trimq)
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("%s: unknown merging metod method" % joining_method)
    qualitycontrol(merged, qc, qcq)
    removechimera(qc, chi, rdb)


def start_otu_pickng(in_folder, out_folder, depth, rdb):
    global PR
    """
    Full analysis will be done
    using bbmerge with maxloose
    quality control with q=19
    using gg database for otu-picking and chimera removal
    set the depth of diversity analysis to 10000
    """

    chi = asfolder(in_folder)
    otus = asfolder(out_folder + PR['Fotus'])
    div = asfolder(out_folder + PR['Fdiv'])

    pickotus(chi, otus, rdb)
    if create_mapping_file:
        create_map(chi, PR['mapping_file'])
    corediv(otus, div, PR['mapping_file'], depth)


def start_diversity_analysis(in_folder, out_folder, mapping_file, depth):
    otus = asfolder(in_folder)
    div = asfolder(out_folder + PR['Fdiv'])
    corediv(in_folder=otus, out_folder=div, mapping_file=mapping_file, depth=depth)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Microbiome analysis using multiple methods
    Version: %s
    Date: %s """ % (__version__, __date__))
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

    parser.add_argument("-b",
                        dest="beginwith",
                        type=str,
                        help="begin with: [otu_picking], [diversity_analysis]")

    parser.add_argument("-t",
                        # "--trim_phred_quality_threshold",
                        dest="trim_threshold",
                        type=int,
                        help="phred quality threshold for trimming [10]",
                        default=10)
    parser.add_argument("-s",
                        dest="stop_at",
                        type=str,
                        help="stop at [chimera_removal\]")

    parser.add_argument("-j",
                        dest='joining_method',
                        help="choose the merging method (fastq-join) or (bbmerge) ",
                        type=str,
                        default="fastq-join")

    parser.add_argument("-d",
                        type=int,
                        dest="fastq_p",
                        help="Percentage of maximum difference in fastq-join [16]",
                        default=16)

    parser.add_argument("-q",
                        dest="qc_threshold",
                        type=int,
                        help="quality control phred threshold [19]",
                        default=19)

    parser.add_argument("-c",
                        dest="ConfigFile",
                        type=str,
                        default='qiime.cfg',
                        help="Configuration file name [qiime.cfg]")

    parser.add_argument("--adapter",
                        dest="adapter_reference",
                        help="Adapters reference file",
                        type=str)
    parser.add_argument("-a",
                        dest="mapping_file",
                        help="Mapping file name",
                        type=str)

    parser.add_argument("-p",
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
                        dest="fungus",
                        help="Using Unite database for fungal samples[False]",
                        action="store_true")

    parser.add_argument("-m",
                        dest="maxloose",
                        help="Assign maxloose to be true for bbmerge [False]",
                        action="store_true")

    parser.add_argument("-r",
                        dest="rdb",
                        help="Reference data base [silva, greengenes]",
                        type=str,
                        default="silva")

    parser.add_argument("-e",
                        dest="depth",
                        type=int,
                        help="set the depth of diversity analyses [10000]",
                        default=10000)

    parser.add_argument("--ml",
                        dest="minimum_length",
                        type=int,
                        help="Minimum length of reads kept after merging [380]",
                        default=380)


    x = parser.format_usage()
    parser.usage = starting_message + x
    arg = parser.parse_args()

    PR.update({
        'in_folder': arg.input,
        'out_folder': arg.output,
        'rdb': arg.rdb,
        'qcq': arg.qc_threshold,
        'maxloose': arg.maxloose,
        'trimq': arg.trim_threshold,
        'joining_method': arg.joining_method,
        'fastq_p': arg.fastq_p,
        'depth': arg.depth,
        'fungus': arg.fungus,
        'ConfigFile': arg.ConfigFile,
        # 'automatic_parameter_file': arg.automatic_parameter_file,
        'parameter_file_name': arg.parameter_file_name,
        # 'reference': arg.reference,
        'beginwith': arg.beginwith,
        'mapping_file': arg.mapping_file,
        'adapter_ref' : arg.adapter_reference,
        'minimum_length': arg.minimum_length})

    ## parameter_file
    get_configuration()
    check_before_start()
    PR['others'] = asfolder(PR['out_folder'] + PR['Fothers'])
    PR['number_of_cores'] = arg.number_of_cores
    ## find if the output folder is present or not
    if (os.path.isdir(PR['out_folder'])):
        sys.exit()
    else:
        os.mkdir(PR['out_folder'])

    if arg.parameter_file_name == None:
        PR['parameter_file_name'] = PR['others'] + "para%s.txt" % PR['id']
        write_parameter_file(PR['parameter_file_name'])
    if arg.mapping_file == None:
        create_mapping_file = True
        PR['mapping_file'] = PR['others'] + "map.tsv"
    else:
        PR['mapping_file'] = arg.mapping_file

    if (arg.beginwith == "diversity_analysis") and (arg.mapping_file == None):
        pass
    if not os.path.isdir(PR['others']):
        os.mkdir(PR['others'])

    logging.basicConfig(filename=PR['others'] + "log.txt",
                        format='%(asctime)s %(levelname)s \n %(message)s',
                        level=logging.DEBUG)

    number_of_cores = PR['number_of_cores']
    loginfo('started')
    [loginfo(str(P) + ": " + str(PR[P])) for P in PR]

    if arg.beginwith == "otu_picking":
        start_otu_pickng(in_folder=PR['in_folder'],
                         out_folder=PR['out_folder'],
                         rdb=PR['rdb'],
                         depth=PR['depth'])

    elif arg.beginwith == "diversity_analysis":
        start_diversity_analysis(in_folder=PR['in_folder'],
                                 out_folder=PR['out_folder'],
                                 depth=PR['depth'])



    elif arg.stop_at == "chimera_removal":
        stop_at_chimera_removal(in_folder=PR['in_folder'],
                                out_folder=PR['out_folder'],
                                rdb=PR['rdb'],
                                joining_method=PR['joining_method'],
                                fastq_p=PR['fastq_p'],
                                maxloose=PR['maxloose'],
                                qcq=PR['qcq'],
                                trimq=PR['trimq'])
    elif arg.stop_at == "merging":
        stop_at_merging(in_folder=PR['in_folder'],
                        out_folder=PR['out_folder'],
                        joining_method=PR['joining_method'],
                        fastq_p=PR['fastq_p'],
                        maxloose=PR['maxloose'],
                        trimq=PR['trimq'])

    elif arg.stop_at == "quality_control":
        stop_at_quality_control(in_folder=PR['in_folder'],
                                out_folder=PR['out_folder'],
                                joining_method=PR['joining_method'],
                                fastq_p=PR['fastq_p'],
                                maxloose=PR['maxloose'],
                                qcq=PR['qcq'],
                                trimq=PR['trimq'])




    else:
        full_analysis(in_folder=PR['in_folder'],
                      out_folder=PR['out_folder'],
                      rdb=PR['rdb'],
                      joining_method=PR['joining_method'],
                      fastq_p=PR['fastq_p'],
                      maxloose=PR['maxloose'],
                      qcq=PR['qcq'],
                      depth=PR['depth'],
                      trimq=PR['trimq'])

    loginfo("Finished")
