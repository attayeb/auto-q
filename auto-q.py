#! /usr/bin/env python
from __future__ import unicode_literals
# coding: utf-8



# Use Gzipped files without extractionpyp


import shutil
import logging
import datetime
import os  # operation system packages
import ConfigParser as configparser  # a package to parse INI file or confige file.
import argparse  # a package to parse commandline arguments.
import sys
from re import sub
import gzip


from subprocess import call  # to run command line scripts
from subprocess import Popen, PIPE, check_output
from multiprocessing.dummy import Pool as Pool

__version__ = '0.2.7.2'
__author__ = "Attayeb Mohsen"
__date__ = "23/1/2019"
## 30/8/2018 add primer-length parameter
## 23/1/2019 correct a merging bug
## add start_at_chimera_removal option

starting_message = """ Microbiome analysis using multiple methods
    Version: %s
    Date: %s
    National Institutes of Biomedical Innovation, Health, and Nutrition\n""" % \
                   (__version__, __date__)

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


def execute(command, shell=True):
    """
    Execute command using os package and return output to log file

    :param command: The command to be executed
    :type command: str
    :param shell: Takes either True or False
    :type shell: boolean
    :return: Run the command in the background and save the
    output to the logging file.

    """
    loginfo(command)
    p = Popen(command.split(), stderr=PIPE, stdout=PIPE)
    output, error = p.communicate()
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
    locate the folder of bbmap
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
                          "please check the configuration file: %s "
                          "to set up the correct folder" % PR['ConfigFile'])
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
    elif PR['rdb'] == "unite":
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
    if os.path.exists(PR['others']):
         pass
    else:
         os.mkdir(PR['others'])

    f = open(parameter_file, "w")
    f.write(parameter_string)
    f.close()

#def copyfilesanddecompress(inFolder, outFolder):
#    shutil.copytree(asfolder(inFolder), asfolder(outFolder))
#    print('copying files')
#    execute("gunzip %s*.gz"%asfolder(outFolder))
#    print('decompress files')

def primertrim(infqfile, outfqfile, length):
    """

    :param infqfile:
    :param outfqfile:
    :param length:
    :return:
    """
    if infqfile.endswith(".gz"):
        infq = gzip.open(infqfile, "r")
    else:
        infq = open(infqfile, "r")
    if outfqfile.endswith(".gz"):
        outfq = gzip.open(outfqfile, "w")
    else:
        outfq = open(outfqfile, "w")
    lines = infq.readlines()
    for a, b, c, d in zip(lines[0::4], lines[1::4], lines[2::4], lines[3::4]):
        outfq.write(a)
        outfq.write(b[length:])
        outfq.write(c)
        outfq.write(d[length:])

    infq.close()
    outfq.close()


def trimfolder(inFolder, outFolder, trimq, ftrim=True):
    """

    """

    import os

    inFolder = asfolder(inFolder)
    outFolder = asfolder(outFolder)

    files = os.listdir(inFolder)
    files.sort()
    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]
    os.mkdir(outFolder)
    # call("mkdir -p %s" % out_folder, shell=True)
    print("Trimming...")

    # get_ipython().system(u'mkdir -p {out_folder}')
    def process(i):
        in1 = inFolder + ins1[i]
        in2 = inFolder + ins2[i]
        print("\n%s and %s" % (ins1[i], ins2[i]))
        out1 = outFolder + ins1[i]
        out2 = outFolder + ins2[i]
        out1_temp1 = outFolder + "temp1_" + ins1[i]
        out2_temp1 = outFolder + "temp1_" + ins2[i]

        # forctrimleft was added
        if ftrim:
            primertrim(in1, out1_temp1, PR['primertrim_forward'])
            primertrim(in2, out2_temp1, PR['primertrim_reverse'])

        else:
            out1_temp1 = in1
            out2_temp1 = in2

        if PR['adapter_ref'] != None:

            execute(
                "bbduk.sh -Xmx1000m -in1=%s -in2=%s -out1=%s -out2=%s -outm=stdout.fa -ref=%s -qtrim=r -trimq=%d -k=18 -ktrim=f" %
                (out1_temp1, out2_temp1, out1, out2, PR['adapter_ref'], trimq), shell=True)
        else:
            execute(
                "bbduk.sh -Xmx1000m -in1=%s -in2=%s -out1=%s -out2=%s -qtrim=r -trimq=%d" %
                (out1_temp1, out2_temp1, out1, out2, trimq), shell=True)

        os.remove(out1_temp1)
        os.remove(out2_temp1)
    p = Pool(PR['number_of_cores'])
    p.map(process, range(len(ins1)))


def mergefolderbb(inFolder, outFolder, maxloose=True):
    """

    """
    inFolder = asfolder(inFolder)
    outFolder = asfolder(outFolder)

    import os
    files = os.listdir(inFolder)
    files.sort()

    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]
    outs = [x.replace("_L001_R1_001", "") for x in ins1]
    os.mkdir(outFolder)
    print("\nMerging ...")

    def process(i):

        in1 = inFolder + ins1[i]
        in2 = inFolder + ins2[i]

        print("%s and %s" % (ins1[i], ins2[i]))
        out = outFolder + outs[i]
        if maxloose:
            execute("bbmerge.sh -in1=%s -in2=%s -out=%s -maxloose=t -ignorebadquality" % (in1, in2, out), shell=True)



        else:
            execute("bbmerge.sh -in1=%s -in2=%s -out=%s -ignorebadquality" % (in1, in2, out), shell=True)

        if PR['remove_intermediate']:
            os.remove(in1)
            os.remove(in2)

    p = Pool(PR['number_of_cores'])
    p.map(process, range(len(ins1)))
    if PR['remove_intermediate']:
        os.removedirs(inFolder)
    print("Merging finished.")


def mergefolder(inFolder, outFolder, pp):
    """

    """
    global PR
    inFolder = asfolder(inFolder)
    outFolder = asfolder(outFolder)

    files = os.listdir(inFolder)

    files.sort()

    ins1 = [x for x in files if "_R1_" in x]
    ins2 = [x.replace("_R1_", "_R2_") for x in ins1]

    outs = [x.replace("_L001_R1_001", "") for x in ins1]
    os.mkdir(outFolder)

    def process(i):
        in1 = inFolder + ins1[i]
        in2 = inFolder + ins2[i]
        print("Merging: %s and %s " % (ins1[i], ins2[i]))
        out = outFolder + "temp_" + outs[i]
        out_final = outFolder + outs[i]
        if out_final.endswith(".gz"):
            out_final = sub(".gz", "", out_final)
        execute("fastq-join -p %d %s %s -o %s" % (pp, in1, in2, out), shell=True)
        os.remove("%sun1" % out)
        os.remove("%sun2" % out)
        os.rename("%sjoin" % out, out)
        remove_short_reads(out, out_final, PR['minimum_length'])
        os.remove(out)
        if PR['remove_intermediate']:
            os.remove(in1)
            os.remove(in2)



    p = Pool(PR['number_of_cores'])
    p.map(process, range(len(ins1)))
    if PR['remove_intermediate']:
        os.removedirs(inFolder)

def qualitycontrol(inFolder, outFolder, q):
    """


    """
    inFolder = asfolder(inFolder)
    outFolder = asfolder(outFolder)

    import os
    files = os.listdir(inFolder)
    files.sort()
    os.mkdir(outFolder)

    # call("mkdir -p %s " % out_folder, shell=True)

    def process(i):
        temp = outFolder + "temp" + i + "/"
        print("\nQuality control: %s" % i)
        sampleId = i.replace(".fastq", "")
        inFile = inFolder + i
        outFile = outFolder + i.replace(".fastq", ".fasta")
        execute("""split_libraries_fastq.py -i %s -o %s --barcode_type not-barcoded --sample_ids %s -q %s""" % (
            inFile, temp, sampleId, q), shell=True)

        tempFile = temp + "seqs.fna"
        call("mv %s %s" % (tempFile, outFile), shell=True)
        call("rm -r %s" % temp, shell=True)
        if PR['remove_intermediate']:
            os.remove(inFile)


    p = Pool(PR['number_of_cores'])
    p.map(process, files)
    print("Quality control finished.")
    if PR['remove_intermediate']:
        os.removedirs(inFolder)

def removechimera(inFolder, outFolder, rdb="silva"):
    """


    """
    global PR
    import os
    inFolder = asfolder(inFolder)
    outFolder = asfolder(outFolder)

    files = os.listdir(inFolder)
    files.sort()

    os.mkdir(outFolder)

    # call("mkdir -p %s" % out_folder, shell=True)

    def process(i):
        print("Chimera removal: %s" % i)
        temp = outFolder + "temp" + i + "/"
        if rdb == "silva":
            execute("identify_chimeric_seqs.py -i %s -m usearch61 -o %s -r %s"
                    % (inFolder + i, temp + i, PR['silva_chim_ref']),
                    shell=True)
        else:
            execute("identify_chimeric_seqs.py -i %s -m usearch61 -o %s -r %s" % (
                inFolder + i, temp + i, PR['gg_chim_ref']),
                    shell=True)

        execute("filter_fasta.py -f %s -o %s -s %s/non_chimeras.txt" % (inFolder + i, outFolder + i, temp + i),
                shell=True)
        call("rm -r %s" % temp, shell=True)
        if PR['remove_intermediate']:
            os.remove(inFolder+i)

    p = Pool(PR['number_of_cores'])
    p.map(process, files)
    if PR['remove_intermediate']:
        os.removedirs(inFolder)

def pickotus(inFolder, outFolder, rdb="silva", fungus=False):
    """


    """

    # TODO : add no parallel option
    global PR
    inFolder = asfolder(inFolder)
    outFolder = asfolder(outFolder)

    inFolder_fasta = inFolder + "*.fasta"
    print("Otu picking...")
    if PR['np']:
        parallel_string = ""
    else:
        parallel_string = "-a -O %d" % PR['number_of_cores']


    if PR['c_ref'] != "none":
        if rdb == "silva":
            execute("pick_open_reference_otus.py -i %s -o %s -p %s -r %s %s -n %s"
                    % (
                        inFolder_fasta, outFolder, PR['parameter_file_name'], PR['c_ref'], parallel_string, PR['c_otu_id']),
                    shell=True)
            #execute("filter_otus_from_otu_table.py -i %s -o %s --negate_ids_to_exclude -e %s"
            #        % (out_folder + "otu_table_mc2_w_tax_no_pynast_failures.biom",
            #           out_folder + "otu_table_mc2_w_tax_no_pynast_failures_close_reference.biom",
            #           PR['silva_reference_seqs']), shell=True)

        elif fungus:
            execute("pick_open_reference_otus.py -i %s -o %s -p %s %s -n %s --suppress_align_and_tree"
                    % (inFolder_fasta, outFolder, PR['parameter_file_name'], parallel_string, PR['c_otu_id']), shell=True)

        else:
            execute("pick_open_reference_otus.py -i %s -o %s -r %s -p %s %s -n %s"
                    % (inFolder_fasta, outFolder,
                       PR['c_ref'], PR['parameter_file_name'],
                       parallel_string, PR['c_otu_id']), shell=True)

            #execute("filter_otus_from_otu_table.py -i %s -o %s --negate_ids_to_exclude -e %s"
            #        % (out_folder + "otu_table_mc2_w_tax_no_pynast_failures.biom",
            #           out_folder + "otu_table_mc2_w_tax_no_pynast_failures_close_reference.biom",
            #           PR['gg_reference_seqs']), shell=True)



    else:
        if rdb == "silva":
            execute("pick_open_reference_otus.py -i %s -o %s -p %s -r %s %s -n %s"
                    % (inFolder_fasta, outFolder, PR['parameter_file_name'], PR['silva_reference_seqs'], parallel_string,
                       PR['c_otu_id']),
                    shell=True)
            execute("filter_otus_from_otu_table.py -i %s -o %s --negate_ids_to_exclude -e %s"
                    % (outFolder + "otu_table_mc2_w_tax_no_pynast_failures.biom",
                       outFolder + "otu_table_mc2_w_tax_no_pynast_failures_close_reference.biom",
                       PR['silva_reference_seqs']), shell=True)

        elif fungus:
            execute("pick_open_reference_otus.py -i %s -o %s -p %s %s -n %s--suppress_align_and_tree"
                    % (inFolder_fasta, outFolder, PR['parameter_file_name'], parallel_string,
                       PR['c_otu_id']), shell=True)

        else:
            execute("pick_open_reference_otus.py -i %s -o %s -r %s -p %s -n %s"
                    % (inFolder_fasta, outFolder,
                       PR['gg_reference_seqs'], PR['parameter_file_name'],
                       parallel_string, PR['c_otu_id']), shell=True)

            execute("filter_otus_from_otu_table.py -i %s -o %s --negate_ids_to_exclude -e %s"
                    % (outFolder + "otu_table_mc2_w_tax_no_pynast_failures.biom",
                       outFolder + "otu_table_mc2_w_tax_no_pynast_failures_close_reference.biom",
                       PR['gg_reference_seqs']), shell=True)

    if PR['remove_intermediate']:
        os.removedirs(inFolder)


def writedf(outFile, ids, sampleIds):
    f = open(outFile, "w+")
    f.write("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tRead\tFile\tDescription\n")
    for x in range(len(ids)):
        f.write("%s\t\t\tR1\t%s\tsingle_file\n" % (ids[x], sampleIds[x]))
    f.close()


def create_map(inFolder, outFile):
    """


    """
    inFolder = asfolder(inFolder)
    print("Writing mapping file")
    import os
    sampleIds = os.listdir(inFolder)
    ids = [x.replace(".fasta", "") for x in sampleIds]
    ids = [x.split("_")[0] for x in ids]
    d = {'#SampleID': ids}
    writedf(outFile, ids, sampleIds)


def corediv(inFolder, outFolder, mappingFile, depth):
    """

    """

    print("Core diversity analyses...")
    inFolder = asfolder(inFolder)
    outFolder = asfolder(outFolder)
    if PR['fungus']:
        biom = inFolder + "otu_table_mc2_w_tax.biom"
    else:
        biom = inFolder + "otu_table_mc2_w_tax_no_pynast_failures.biom"
    tree = inFolder + "rep_set.tre"
    # get_ipython().system(
    #    u'core_diversity_analyses.py -i {biom}     -o {out_folder}     -m {mapping_file}     -t {tree}     -e {depth}')
    if PR['fungus']:
        execute("core_diversity_analyses.py -i %s -o %s -m %s -e %d --nonphylogenetic_diversity" % (
            biom, outFolder, mappingFile, depth),
                shell=True)
    else:
        execute(
            "core_diversity_analyses.py -i %s -o %s -m %s -t %s -e %d" % (biom, outFolder, mappingFile, tree, depth),
            shell=True)


def full_analysis(inFolder, outFolder, depth, rdb, trimq, joining_method,
                  qcq, maxloose, fastq_p):
    global PR
    """

    """
    trimmed = asfolder(outFolder + PR['Ftrimmed'])
    merged = asfolder(outFolder + PR['Fmerged'])
    qc = asfolder(outFolder + PR['Fqc'])
    chi = asfolder(outFolder + PR['Fchi'])
    otus = asfolder(outFolder + PR['Fotus'])
    div = asfolder(outFolder + PR['Fdiv'])

    trimfolder(inFolder, trimmed, trimq)
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("Wrong method")
    qualitycontrol(merged, qc, qcq)
    removechimera(qc, chi, rdb)
    pickotus(chi, otus, rdb)
    # here
    if create_mapping_file:
        create_map(qc, PR['mapping_file'])
    corediv(otus, div, PR['mapping_file'], depth)


def stop_at_merging(inFolder, outFolder, trimq, joining_method, maxloose, fastq_p):
    global PR
    trimmed = asfolder(outFolder + PR['Ftrimmed'])
    merged = asfolder(outFolder) + PR['Fmerged']
    trimfolder(inFolder, trimmed, trimq)
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("%s: unknown merging metod method" % joining_method)


def stop_at_quality_control(inFolder, outFolder, joining_method, trimq,
                            qcq, maxloose, fastq_p):
    global PR
    """
    """

    trimmed = asfolder(outFolder + PR['Ftrimmed'])
    merged = asfolder(outFolder + PR['Fmerged'])
    qc = asfolder(outFolder + PR['Fqc'])

    trimfolder(inFolder, trimmed, trimq)
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("%s: unknown merging metod method" % joining_method)
    qualitycontrol(merged, qc, qcq)


def stop_at_chimera_removal(inFolder, outFolder, rdb, trimq, joining_method,
                            qcq, maxloose, fastq_p):
    """

    """

    global PR

    trimmed = asfolder(outFolder + PR['Ftrimmed'])
    merged = asfolder(outFolder + PR['Fmerged'])
    qc = asfolder(outFolder + PR['Fqc'])
    chi = asfolder(outFolder + PR['Fchi'])

    trimfolder(inFolder, trimmed, trimq)
    if joining_method == "fastq-join":
        mergefolderfastq(trimmed, merged, fastq_p)
    elif joining_method == "bbmerge":
        mergefolderbb(trimmed, merged, maxloose=maxloose)
    else:
        raise ("%s: unknown merging metod method" % joining_method)
    qualitycontrol(merged, qc, qcq)
    removechimera(qc, chi, rdb)


def start_at_chimera_removal(inFolder, outFolder, rdb, depth):
    global PR

    qc = asfolder(inFolder)
    chi = asfolder(outFolder + PR['Fchi'])
    otus = asfolder(outFolder + PR['Fotus'])
    div = asfolder(outFolder + PR['Fdiv'])

    removechimera(qc, chi, rdb)
    pickotus(chi, otus, rdb)
    # here
    if create_mapping_file:
        create_map(qc, PR['mapping_file'])
    corediv(otus, div, PR['mapping_file'], depth)

def start_otu_pickng(inFolder, outFolder, depth, rdb):
    """

    """
    global PR
    chi = asfolder(inFolder)
    otus = asfolder(outFolder + PR['Fotus'])
    div = asfolder(outFolder + PR['Fdiv'])

    pickotus(chi, otus, rdb)

    if create_mapping_file:
        create_map(chi, PR['mapping_file'])

    corediv(otus, div, PR['mapping_file'], depth)


def start_diversity_analysis(inFolder, outFolder, mapping_file, depth):
    otus = asfolder(inFolder)
    div = asfolder(outFolder + PR['Fdiv'])
    corediv(inFolder=otus, outFolder=div, mappingFile=mapping_file, depth=depth)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Microbiome analysis using multiple methods
    Version: %s
    Date: %s """ % (__version__, __date__))
    parser.add_argument("-i",  # "--input",
                        dest="input",
                        # type=str,
                        help="the input sequences filepath (fastq files) [REQUIRED]",
                        metavar="Input folder",
                        required=True)

    parser.add_argument("-o",
                        # "--output",
                        dest="output",
                        type=str,
                        metavar="Output folder",
                        help="the output directory [REQUIRED]",
                        required=True)


    parser.add_argument("-t",
                        dest="trim_threshold",
                        type=int,
                        metavar="trim_phred_threshold",
                        help="phred quality threshold for trimming [default: 12]",
                        default=12)

    parser.add_argument("-p",
                        type=int,
                        dest="fastq_p",
                        metavar="fastq-join p",
                        help="fastq-join's percentage of mismatch [default: 16]",
                        default=16)

    parser.add_argument("--adapter",
                        metavar=None,
                        dest="adapter_reference",
                        help="Adapters reference file",
                        type=str)

    parser.add_argument("-b",
                        dest="beginwith",
                        type=str,
                        metavar="starting step",
                        choices=['otu_picking', 'diversity_analysis', 'chimera_removal'],
                        help="starting the analysis in the middle: (otu_picking), (diversity_analysis), (chimera_removal)")


    parser.add_argument("-s",
                        dest="stop_at",
                        type=str,
                        metavar="stop at",
                        choices = ['merging', 'quality_control','chimera_removal'],
                        help='terminate the analysis at this step [choices: (merging), (quality_control), (chimera_'
                             'removal))')

    parser.add_argument("-j",
                        dest='joining_method',
                        help="choose the merging method (fastq-join) or (bbmerge) [default: fastq-join]",
                        type=str,
                        metavar="joining method",
                        choices = ['fastq-join', "bbmerge"],
                        default="fastq-join")

    parser.add_argument("-m",
                        dest="maxloose",
                        help="Assign maxloose to be true for bbmerge [default: False]",
                        action="store_true")

    parser.add_argument("-q",
                        dest="qc_threshold",
                        type=int,
                        metavar="quality control threshold",
                        help="quality control phred threshold [default: 19]",
                        default=19)

    parser.add_argument("--continuation_reference",
                        dest="c_ref",
                        type=str,
                        metavar="newref_seq.fna",
                        help="reference sequence for continuation. If you want to continue analysis using the reference "
                             "data set from previous analysis. you can find it in the last sample otus folder new_refseqs.fna",
                        default="none")

    parser.add_argument("--continuation_otu_id",
                        dest="c_otu_id",
                        type=str,
                        metavar=None,
                        help="continuation reference new otus ids",
                        default="New")

    parser.add_argument("-r",
                        dest="rdb",
                        metavar="Reference database",
                        help="silva, greengenes [default: silva]",
                        choices=['silva', 'greengenes', 'unite'],
                        type=str,
                        default="silva")

    parser.add_argument("-c",
                        dest="ConfigFile",
                        type=str,
                        metavar="Configuration file name",
                        default='qiime.cfg',
                        help="Configuration file name [default: qiime.cfg]")

    parser.add_argument("-a",
                        dest="mapping_file",
                        help="Mapping file name",
                        metavar="Mapping file name",
                        type=str)

    parser.add_argument("--parameter_file_name",
                        help="The name of the parameter file [if not assigned is automatically produced using "
                             "configuration file",
                        type=str,
                        metavar=None,
                        dest="parameter_file_name")

    parser.add_argument("-n",
                        # "--number_of_cores",
                        help="Specify the number of jobs to start with [default: 2]",
                        type=int,
                        metavar='Number of jobs',
                        dest="number_of_cores",
                        default=2)

    parser.add_argument("-e",
                        dest="depth",
                        type=int,
                        metavar="Sampling depth",
                        help="sampling depth for diversity analyses [default: 10000]",
                        default=10000)

    parser.add_argument("--remove_intermediate_files",
                        help="To remove intermediate files, to reduce the disk space",
                        dest="remove_intermediate",
                        action="store_true")

#    parser.add_argument("--decompress",
#                        help="Copy input files to outputfolder/fastq and decompress them",
#                        dest="decompress",
#                        action="store_true")

    parser.add_argument("--ml",
                        dest="minimum_length",
                        metavar='Minimum length',
                        type=int,
                        help="Minimum length of reads kept after merging [default: 380]",
                        default=380)

    parser.add_argument("--primer-trim-f",
                        dest="primertrim_forward",
                        metavar='Primer Trim',
                        type=int,
                        help="length of the forward primer [17]",
                        default=17)

    parser.add_argument("--primer-trim-r",
                        dest="primertrim_reverse",
                        metavar='Primer Trim',
                        type=int,
                        help="length of the reverse primer [21]",
                        default=21)

    #x = parser.format_usage()
    #parser.usage = starting_message #+ x
    arg = parser.parse_args()

    PR.update({

        'in_folder': asfolder(arg.input),
        'out_folder': asfolder(arg.output),
#        'decompress': arg.aaa
        # ress,
        'rdb': arg.rdb,
        'qcq': arg.qc_threshold,
        'maxloose': arg.maxloose,
        'trimq': arg.trim_threshold,
        'joining_method': arg.joining_method,
        'fastq_p': arg.fastq_p,
        'depth': arg.depth,
        'ConfigFile': arg.ConfigFile,
        'parameter_file_name': arg.parameter_file_name,
        'remove_intermediate': arg.remove_intermediate,
        'beginwith': arg.beginwith,
        'mapping_file': arg.mapping_file,
        'adapter_ref': arg.adapter_reference,
        'minimum_length': arg.minimum_length,
        'c_ref': arg.c_ref,
        'c_otu_id': arg.c_otu_id,
        'primertrim_forward': arg.primertrim_forward,
        'primertrim_reverse': arg.primertrim_reverse})

    ## parameter_file
    get_configuration()
    check_before_start()






    if PR['rdb'] == 'unite':
        PR['fungus'] = True
    else:
        PR['fungus'] = False
    PR['others'] = asfolder(PR['out_folder'] + PR['Fothers'])
    PR['number_of_cores'] = arg.number_of_cores
    if PR['number_of_cores'] == 1:
        PR['np'] = True
    else:
        PR['np'] = False

    if (os.path.isdir(PR['out_folder'])):
        sys.exit()
    else:
        os.mkdir(PR['out_folder'])
    if not os.path.isdir(PR['others']):
        os.mkdir(PR['others'])

    logging.basicConfig(filename=PR['others'] + "log.txt",
                        format='%(levelname)s \n %(message)s',
                        level=logging.DEBUG)
    loginfo('started')
    [loginfo(str(P) + ": " + str(PR[P])) for P in PR]

#    if PR['decompress']:
#        copyfilesanddecompress(PR['in_folder'], asfolder(PR['out_folder']+"fastq"))
#        PR['in_folder'] = asfolder(PR['out_folder'])+'fastq/'


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

    number_of_cores = PR['number_of_cores']


    if arg.beginwith == "otu_picking":
        start_otu_pickng(inFolder=PR['in_folder'],
                         outFolder=PR['out_folder'],
                         rdb=PR['rdb'],
                         depth=PR['depth'])

    elif arg.beginwith == "diversity_analysis":
        start_diversity_analysis(inFolder=PR['in_folder'],
                                 outFolder=PR['out_folder'],
                                 mapping_file=PR['mapping_file'],
                                 depth=PR['depth'])

    elif arg.beginwith == "chimera_removal":
        start_at_chimera_removal(inFolder=PR['in_folder'],
                                 outFolder=PR['out_folder'],
                                 rdb= PR['rdb'],
                                 depth=PR['depth'])

    elif arg.stop_at == "chimera_removal":
        stop_at_chimera_removal(inFolder=PR['in_folder'],
                                outFolder=PR['out_folder'],
                                rdb=PR['rdb'],
                                joining_method=PR['joining_method'],
                                fastq_p=PR['fastq_p'],
                                maxloose=PR['maxloose'],
                                qcq=PR['qcq'],
                                trimq=PR['trimq'])
    elif arg.stop_at == "merging":
        stop_at_merging(inFolder=PR['in_folder'],
                        outFolder=PR['out_folder'],
                        joining_method=PR['joining_method'],
                        fastq_p=PR['fastq_p'],
                        maxloose=PR['maxloose'],
                        trimq=PR['trimq'])

    elif arg.stop_at == "quality_control":
        stop_at_quality_control(inFolder=PR['in_folder'],
                                outFolder=PR['out_folder'],
                                joining_method=PR['joining_method'],
                                fastq_p=PR['fastq_p'],
                                maxloose=PR['maxloose'],
                                qcq=PR['qcq'],
                                trimq=PR['trimq'])




    else:
        full_analysis(inFolder=PR['in_folder'],
                      outFolder=PR['out_folder'],
                      rdb=PR['rdb'],
                      joining_method=PR['joining_method'],
                      fastq_p=PR['fastq_p'],
                      maxloose=PR['maxloose'],
                      qcq=PR['qcq'],
                      depth=PR['depth'],
                      trimq=PR['trimq'])

    loginfo("Finished")

