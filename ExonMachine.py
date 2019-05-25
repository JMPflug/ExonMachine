#!/usr/bin/env python2.7

from __future__ import division
from contextlib import contextmanager
import sys
import re
import numpy as np
from Bio import SeqIO
from Bio import SearchIO
from Bio import AlignIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
import multiprocessing as mp
import subprocess
import os
import argparse
import time
import shlex
import shutil

parser = argparse.ArgumentParser(description="Generates a list of candidate probe loci from a set of aligned or unaligned\
                                            transcriptome sequences in a specified folder using a reference genome to\
                                            determine exon bondaries. Transcriptome sequences should be provided as\
                                            individual multi-fasta files, one file per ortholog, in a folder in the current\
                                            working directory.\
                                            Dependencies: Exonerate, CD-HIT-EST Biopython")\


parser.add_argument("-reference_taxon", "-t", required=True, help="Specify which taxon in your ortholog clusters to use as the source for reference sequences.",
                    type=str)

parser.add_argument("-reference_genome", "-g", required=True, help="Specify a fasta file to use as the reference genome.",
                    type=str)

parser.add_argument("-ortho_directory", "-o", required=True, help="Directory containing ortholog fasta files.",
                    type=str)

parser.add_argument("-T", "--threads", default=1, help="Number of threads to use for exonerate multithreading (default=1).",
                    type=int)


##### Paths to executables (if not already in PATH) ####
transdecoder_LO_path = "/Volumes/Basement_Shelves/EverythingisTerrible/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs"
transdecoder_Pd_path = "/Volumes/Basement_Shelves/EverythingisTerrible/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict"
path_macse = "/Applications/macse_v1.2.jar"
path_mafft = "/usr/local/bin/mafft"
#####################################################################

args = parser.parse_args()
reftrantaxa = args.reference_taxon
refgenome = args.reference_genome
threads = args.threads
orthofolder = args.ortho_directory
cwd = os.getcwd()

def slistdir(directory):
    """A specialized version of os.listdir() that ignores files that
    start with a leading period."""
    filelist = os.listdir(directory)
    return [x for x in filelist
        if not (x.startswith('.'))]

@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def clean_up_working(exon_folder="None"):
    """Deletes all files and directories in current directory except for python scripts. Excludes folder of Exonerate alignments, if provided"""
    for files in slistdir("./"):
        if exon_folder != "None":
            if ".py" not in files and exon_folder.replace(".","").replace("/","") not in files:
                if os.path.isfile(files) == True:
                    os.remove(files)
                if os.path.isdir(files):
                    shutil.rmtree(files)
        else:
            if ".py" not in files:
                if os.path.isfile(files) == True:
                    os.remove(files)
                if os.path.isdir(files):
                    shutil.rmtree(files)

def first_exonerate_parse(dir, newdir, prefix):
    print "first_exonerate_parse"
    cwd = os.getcwd()
    if not os.path.exists(cwd + newdir):
        os.makedirs(cwd + newdir)
    for file in slistdir(cwd + dir):
        if 'DS_Store' not in file:
            result = SearchIO.parse(cwd + dir + file, 'exonerate-text')
            for h in result:
                for hh in h:
                    for hhh in hh:
                        hitcounter = 1
                        for hhhh in hhh:
                            hitseq =  hhhh.query
                            rootname = file.split('.fasta')
                            orthoname = file.split("_")
                            orthosubdir = cwd + newdir + '/' + orthoname[0]
                            if not os.path.exists(orthosubdir):
                                os.makedirs(orthosubdir)
                            newseqstr = str(hitseq.seq.ungap("-"))
                            newid = prefix + str(hitcounter) + '_' + rootname[0]
                            record = SeqRecord(Seq(newseqstr, generic_dna), id =  newid, description = '')
                            fastaname = prefix + str(hitcounter) + '_' + rootname[0] + '.fasta'
                            SeqIO.write(record, orthosubdir + '/' + fastaname, "fasta")
                            hitcounter += 1
    #sys.exit("Done")

### CLEAN DIRECTORY REFERENCES!!! ###
def rename_exons():
    print "rename_exons"
    for orthofolder in os.listdir('./ref_exons/'):
        if "DS_Store" not in orthofolder:
            OFloc = './ref_exons/' + orthofolder
            for ortho in os.listdir(OFloc):
                if '.fa' in ortho and 'Locus' not in ortho:
                    with open('locus_log.txt') as f:
                        for line in f:
                            locuslist = line.split(' ')
                            seqname = re.split('Ortho\d+_', ortho)
                            seqfind = seqname[1].replace('.fasta', '')
                            lls = locuslist[1].strip()
                            if lls == seqfind:
                                newname = OFloc + '/' + orthofolder + '_' + seqname[0].replace('Ref','') + 'Reference.fasta'
                                oldname = OFloc + '/' + ortho
                                os.rename(oldname, newname)

### CLEAN DIRECTORY REFERENCES!!! ###
def ungap_split(in_dir):
    print "ungap_split"
    for file in os.listdir(in_dir):
        if "Ortho" in file:
            for locusfolder in os.listdir('./split_mafft_fasta/' + file):
                for seq_record in SeqIO.parse('./split_mafft_fasta/' + file + '/' + locusfolder, "fasta"):
                    newseq = seq_record.seq.ungap("-")
                    os.remove('./split_mafft_fasta/' + file + '/' + locusfolder)
                    with open('./split_mafft_fasta/' + file + '/' + locusfolder, "w+") as f:
                        f.write('>' + locusfolder + '\n')
                        f.write(str(newseq.upper()))

### CLEAN DIRECTORY REFERENCES!!! ###
def multirun_exonerate(refgenome):
    print "multirun_exonerate"
    if not os.path.exists('./t2g_exn_files/'):
        os.makedirs('./t2g_exn_files/')
    commands = []
    for ortho_file in os.listdir('./ref_orthos/'):
        if 'DS_Store' not in ortho_file:
            ffile = './ref_orthos/' + ortho_file
            outfile = './t2g_exn_files/' + ortho_file + '.txt'
            outfile_name = ortho_file + '.txt'
            #command = 'exonerate ' + "'" + ffile + "' " + refgenome + ' -n ' + '1 ' + '-m ' + 'est2genome ' + '> ' + "'" + outfile + "'"
            command = 'exonerate ' + "'" + ffile + "' " + refgenome + ' -n ' + '1 ' + '--ryo "PercentId_%%ps_PercentSim_%%pi_Query_%%qi \n" ' + '> ' + "'" + outfile + "'"

            if outfile_name not in slistdir('./t2g_exn_files/'):
                print "Running exonerate for " + str(ortho_file)
                commands.append(command)
            else:
                print str(outfile) + " exists!"
    return commands

########## Run Exonerate using "querychunk" function ############################################################

def append_ortho_id_to_fastas(in_dir, out_file_name):
    print "merge_fastas"
    locus_log = locus_log_to_dict("locus_log.txt")
    ortho_count = 1
    if os.path.exists(out_file_name):
        os.remove(out_file_name)
    for file_name in slistdir(in_dir):
        file_name = os.path.join(in_dir, file_name)
        print file_name
        for seq_record in SeqIO.parse(file_name, "fasta"):
            seq_record.id =    "%s_%s" % (locus_log[seq_record.id], str(seq_record.id))
            ortho_count = ortho_count + 1
            seq_record.description = ""
            with open(out_file_name, "a+") as f:
                SeqIO.write(seq_record, f, "fasta")

def server_run_exonerate(in_dir, reference, out_dir, number_results):
    print "server_run_exonerate"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    commands = []
    ortho_count = len(slistdir(in_dir))
    for ortholog in slistdir(in_dir):
        ortho_number = ortholog.split("_")[0]
        ortho_full_path = os.path.join(in_dir, ortholog)
        outfile_name = "Exonerate_outfile_%s.txt" % (ortho_number)
        outfile_full = os.path.join(out_dir, outfile_name)
        outfile_check = "Exonerate_outfile_%s_clean.txt" % (ortho_number)
        outfile_check_full = os.path.join(out_dir, outfile_check)
        if outfile_name not in slistdir(out_dir) and outfile_check not in slistdir(out_dir):
            command = "exonerate  %s %s -m est2genome -n %s  > %s" % (ortho_full_path, reference, number_results, outfile_full)
            print "Preparing to run Exonerate on %s" % (ortho_number)
            commands.append(command)
        else:
            print "Exonerate file for %s already exists! Skipping..." % (ortho_number)
    print  commands
    return commands

def server_run_exonerate_paralogs(in_dir, reference, out_dir, number_results):
    print "server_run_exonerate"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    commands = []
    ortho_count = len(slistdir(in_dir))
    for ortholog in slistdir(in_dir):
        ortho_number = ortholog.split("_")[0]
        exon_number = ortholog.split("_")[1]
        ortho_full_path = os.path.join(in_dir, ortholog)
        outfile_name = "Exonerate_outfile_%s_%s.txt" % (ortho_number, exon_number)
        outfile_full = os.path.join(out_dir, outfile_name)
        outfile_check = "Exonerate_outfile_%s_%s_clean.txt" % (ortho_number, exon_number)
        outfile_check_full = os.path.join(out_dir, outfile_check)
        if outfile_name not in slistdir(out_dir) and outfile_check not in slistdir(out_dir):
            command = "exonerate  %s %s -m affine:local -n %s  > %s" % (ortho_full_path, reference, number_results, outfile_full)
            print "Preparing to run Exonerate on %s" % (ortho_number)
            commands.append(command)
        else:
            print "Exonerate file for %s already exists! Skipping..." % (ortho_number)
    print  commands
    pool = mp.Pool(int(threads))
    pool.map(multihold, commands)

def parallel_server_exonerate(in_dir, reference, out_dir, number_results):
    print "parallel_server_exonerate"
    commands = server_run_exonerate(in_dir, reference, out_dir, number_results)
    pool = mp.Pool(int(threads))
    pool.map(multihold, commands)

def multirun_find_in_alignment(ref_folder, nonref_folder, probe_out_folder):
    print "multirun_find_in_alignment"
    if not os.path.exists(probe_out_folder):
        os.makedirs(probe_out_folder)
    commands = []
    for ortho_file in slistdir(nonref_folder):
        print "\t" + ortho_file
        if 'Ortho' in ortho_file:
            for nonref_ortho in slistdir(nonref_folder + ortho_file):
                for ref_ortho in slistdir(ref_folder):
                    ref_ortho_num = ref_ortho.split('_')[0]
                    nonref_ortho_num = nonref_ortho.split('_')[0]
                    if ref_ortho_num == nonref_ortho_num:
                        query_input = os.path.join(ref_folder, ref_ortho)
                        target_input = os.path.join(nonref_folder, ortho_file, nonref_ortho)
                        taxon_name = re.sub(r'Ortho[0-9]+_', '', nonref_ortho)
                        exonerate_out_name = str(ref_ortho.replace("_Reference.fasta",""))
                        probe_out_path = os.path.join(probe_out_folder, exonerate_out_name)
                        if os.path.exists(probe_out_path):
                            pass
                        else:
                            os.makedirs(probe_out_path)
                        exonerate_out = os.path.join(probe_out_folder, exonerate_out_name, taxon_name + ".txt")
                        #command = 'exonerate ' + " '" + query_input + "' '" + target_input + "'" + ' --ryo "\%%qi_PercentID:_%%ps"\ -n 1 -m coding2coding > ' + "'" + exonerate_out + "'"
                        #command = 'exonerate  %s %s -n 1 -m coding2coding --ryo "PercentId_%%ps_PercentSim_%%pi_Query_%%qi \n"  > %s' % (query_input, target_input, exonerate_out)
                        command = 'exonerate  %s %s -n 1 --ryo "PercentId_%%ps_PercentSim_%%pi_Query_%%qi \n"  > %s' % (query_input, target_input, exonerate_out)
                        #command = 'exonerate ' + " '" + drffile + "' '" + drsubfile + "'" + ' -n 1 -m coding2coding --ryo "Pid%pi_Psm%ps" > ' + "'" + outfile + "'"
                        commands.append(command)
    return commands

def parallel_split_exonerate(in_fasta):
    print "paralle_split_exonerate"
    commands = split_run_exonerate(refgenome, in_fasta)
    pool = mp.Pool(int(threads))
    pool.map(multihold, commands)

#split_exonerate_parse('/t2g_exn_files/', '/ref_exons/', 'Exon')
def split_exonerate_parse(in_dir, new_dir, prefix, sim_threshold):
    clean_up_exonerate_files(in_dir)
    print "split_exonerate_parse"
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    else:
        shutil.rmtree(new_dir)
        os.makedirs(new_dir)
    for file_name in slistdir(in_dir):
        exonerate_file = os.path.join(in_dir,file_name)
        result = SearchIO.parse(exonerate_file, 'exonerate-text')
        print_ortho_switch = True
        for h in result:
            hitcounter = 0
            for hh in h:
                for hhh in hh:
                    for hhhh in hhh:
                        hitcounter += 1
                        hitseq =  hhhh.query
                        sim_string = hhhh.aln_annotation['similarity']
                        len_sim_string = len(sim_string)
                        match_count = sim_string.count('|')
                        sim_score = float(match_count/len_sim_string)
                        if print_ortho_switch == True:
                            print hhhh.query.id.split("_")[0]
                            print_ortho_switch = False
                        if sim_score >= sim_threshold:
                            orthoname = hitseq.id.split("_")[0]
                            new_ortho_name = "%s_%s%s_Reference" % (orthoname, prefix, str(hitcounter))
                            new_file_name =  "%s.fasta" % (new_ortho_name)
                            hitseq.id = new_ortho_name
                            hitseq.seq = hitseq.seq.ungap("-")
                            out_file_path = os.path.join(new_dir, new_file_name)
                            SeqIO.write(hitseq, out_file_path, "fasta")
                            print "\t%s passed (Sim: %s)" % (hhhh.query.id, sim_score)
                        else:
                            orthoname = hitseq.id.split("_")[0]
                            new_ortho_name = "%s_%s%s_Reference" % (orthoname, prefix, str(hitcounter))
                            hitseq.id = new_ortho_name
                            print "      X %s too divergent (Sim: %s)" % (hhhh.query.id, sim_score)

def remove_paralogs(in_dir, exon_seq_dir, prefix, sim_threshold, len_threshold):
    clean_up_exonerate_files(in_dir)
    print "remove_paralogs"
    def find_paralogs():
        exonerate_files_for_removal = []
        for file_name in slistdir(in_dir):
            exonerate_file = os.path.join(in_dir,file_name)
            result = SearchIO.parse(exonerate_file, 'exonerate-text')
            print_ortho_switch = True
            for h in result:
                ortho_fail_switch = False
                sim_threshold_counter = 0
                score_list = []
                #print h
                hitcounter = 0
                query_id = h[0][0].query.id
                print query_id
                hit_list = h
                for hit in hit_list:
                    entry_dict = {}
                    hitcounter += 1
                    match = hit[0]
                    sim_string = match.aln_annotation['similarity']
                    len_sim_string = len(sim_string)
                    match_count = sim_string.count('|')
                    sim_score = float(match_count/len_sim_string) #
                    len_score = len(match.hit)                      #
                    if sim_score >= sim_threshold:
                        sim_threshold_counter += 1
                    if sim_threshold_counter >= 2 and sim_score >= sim_threshold and len_score >= len_threshold:
                        print "      X %s: Sim: %s Len: %s Off-target hit detected!" % (hitcounter, sim_score, len_score)
                        ortho_fail_switch = True
                    else:
                        print "\t%s: Sim: %s Len: %s" % (hitcounter, sim_score, len_score)
                if ortho_fail_switch == False:
                    print "\t%s passed" % (query_id)
                else:
                    print "\t%s failed; %s removed from further consideration" % (query_id, exonerate_file)
                    exonerate_files_for_removal.append(query_id)
        return exonerate_files_for_removal

    def move_failed_exonerate_files(in_list, new_dir):
        if not os.path.exists(new_dir):
            os.makedirs(new_dir)
        else:
            shutil.rmtree(new_dir)
            os.makedirs(new_dir)
        for failed_file in in_list:
            failed_fasta = failed_file + ".fasta"
            failed_file_path = os.path.join(exon_seq_dir, failed_fasta)
            print "Moving %s" % (failed_file_path)
            #file_name = failed_file.split("/")[2]
            new_path = os.path.join(new_dir, failed_fasta)
            os.rename(failed_file_path, new_path)

    move_failed_exonerate_files(find_paralogs(), "./paralogs/")

                #print hit
                #print "\n"

            # hhhh =  h[0][0][0]
            # print hhhh
            # hitcounter += 1
            # hitseq =  hhhh.query
            # sim_string = hhhh.aln_annotation['similarity']
            # len_sim_string = len(sim_string)
            # match_count = sim_string.count('|')
            # sim_score = float(match_count/len_sim_string)
            # if print_ortho_switch == True:
            #     print hhhh.query.id.split("_")[0]
            #     print_ortho_switch = False



            # if sim_score >= sim_threshold:
            #     orthoname = hitseq.id.split("_")[0]
            #     new_ortho_name = "%s_%s%s_Reference" % (orthoname, prefix, str(hitcounter))
            #     new_file_name =  "%s.fasta" % (new_ortho_name)
            #     hitseq.id = new_ortho_name
            #     out_file_path = os.path.join(new_dir, new_file_name)
            #     SeqIO.write(hitseq, out_file_path, "fasta")
            #     print "\t%s passed (Sim: %s)" % (hhhh.query.id, sim_score)
            # else:
            #     orthoname = hitseq.id.split("_")[0]
            #     new_ortho_name = "%s_%s%s_Reference" % (orthoname, prefix, str(hitcounter))
            #     hitseq.id = new_ortho_name
            #     print "      X %s too divergent (Sim: %s)" % (hhhh.query.id, sim_score)

def locus_log_to_dict(log_file):
    locus_log_dict = {}
    with open(log_file, "r+") as f:
        for line in f:
            line_split = line.strip().split()
            key = line_split[1]
            if key not in locus_log_dict:
                locus_log_dict[key] = line_split[0]
    return locus_log_dict

def merge_fastas(in_dir, out_file_name):
    print "merge_fastas"
    locus_log = locus_log_to_dict("locus_log.txt")
    ortho_count = 1
    if os.path.exists(out_file_name):
        os.remove(out_file_name)
    for file_name in slistdir(in_dir):
        file_name = os.path.join(in_dir, file_name)
        print file_name
        for seq_record in SeqIO.parse(file_name, "fasta"):
            seq_record.id =    "%s_%s" % (locus_log[seq_record.id], str(seq_record.id))
            ortho_count = ortho_count + 1
            seq_record.description = ""
            with open(out_file_name, "a+") as f:
                SeqIO.write(seq_record, f, "fasta")

def clean_up_exonerate_files(in_dir):
    print "clean_up_exonerate_files"
    for file_name in slistdir(in_dir):
        if "_clean.txt" not in file_name:
            in_file = os.path.join(in_dir, file_name)
            with open(in_file, "r+") as f:
                for line in f:
                    new_line = line.replace("[revcomp]","")
                    new_file_name = in_file.replace(".txt", "_clean.txt")
                    with open(new_file_name, "a+") as g:
                        g.write(new_line)
            os.remove(in_file)

def qc_exon_to_fasta(in_dir, out_dir, sim_threshold):
    print "exon_to_fasta"
    exonerate_probes_folder = in_dir
    probe_seq_folder = out_dir
    if not os.path.exists(probe_seq_folder):
        os.makedirs(probe_seq_folder)
    for root, dirs, exonerate_files in os.walk(exonerate_probes_folder, topdown=False):
        print_ortho_switch = True
        exonerate_files = [f for f in exonerate_files if not f[0] == '.']
        for exonerate_file in exonerate_files:
            exonerate_file_path = os.path.join(root, exonerate_file)
            result = SearchIO.parse(exonerate_file_path, 'exonerate-text')
            for h in result:
                hit =  h[0][0][0]
                hitseq =  hit.hit
                ortho_len = len(hitseq.seq)
                hitseq.id = hitseq.id.replace('.fasta', '')
                sim_string = hit.aln_annotation['similarity']
                len_sim_string = len(sim_string)
                match_count = sim_string.count('|')
                sim_score = float(match_count/len_sim_string)
                hitseq.description = "len: " + str(ortho_len)
                hitseq.seq = hitseq.seq.ungap("-")
                ortho_name = hit.query.id.replace("_Reference","")
                taxon_name = re.sub(r'Ortho[0-9]+_', '', hitseq.id)
                full_exon_name = ortho_name + "_" + taxon_name
                hitseq.id = full_exon_name
                fastaname = "%s_len_%s_%s.fasta" % (ortho_name, str(ortho_len), taxon_name)
                probe_seq_dir = os.path.join(probe_seq_folder, ortho_name)
                if not os.path.exists(probe_seq_dir):
                    os.makedirs(probe_seq_dir)
                probe_seq_out_path = os.path.join(probe_seq_dir, full_exon_name + ".fasta")
                if print_ortho_switch == True:
                    print ortho_name
                    print_ortho_switch = False
                if sim_score >= sim_threshold:
                    SeqIO.write(hitseq, probe_seq_out_path, "fasta")
                    print "\t%s passed (Sim: %s)" % (hitseq.id, sim_score)
                else:
                    print "      X %s too divergent (Sim: %s)" % (hitseq.id, sim_score)

### Probably okay to delete ### CLEAN DIRECTORY REFERENCES!!! ###
def exon_to_fasta():
    print "exon_to_fasta"
    cwd = os.getcwd()
    if not os.path.exists('./probe_seqs'):
        os.makedirs('./probe_seqs')
    if not os.path.exists('./probes_merged/'):
        os.makedirs('./probes_merged/')
    for ortho_file in slistdir('./probes'):
        if 'Ortho' in ortho_file:
            print "\t" + ortho_file
            seq_list = []
            for subfile in os.listdir('./probes/' + ortho_file):
                if 'Ortho' in subfile:
                    EL = './probes/' + ortho_file + '/'
                    os.chdir(cwd + '/probes/' + ortho_file)
                    result = SearchIO.parse(subfile, 'exonerate-text')
                    for h in result:
                        hit =  h[0][0][0]
                        hitseq =  hit.hit
                        ortho_len = len(hitseq.seq)
                        hitseq.id = hitseq.id.replace('.fasta', '')
                        hitseq.description = "len: " + str(ortho_len)
                        taxon_name = re.sub(r'Ortho[0-9]+_', '', subfile.replace('.fasta',''))
                        taxon_name = taxon_name.split('_vs_')[1]
                        fastaname = ortho_file + "_len_" + str(ortho_len) + '_' + taxon_name
                        hitseq.seq = hitseq.seq.ungap("-")
                        if not os.path.exists(cwd + '/probe_seqs/' + ortho_file):
                            os.makedirs(cwd + '/probe_seqs/' + ortho_file)
                        SeqIO.write(hitseq, cwd + '/probe_seqs/' + ortho_file + '/' + fastaname, "fasta")
                        seq_list.append(hitseq)
                    os.chdir(cwd)
            SeqIO.write(seq_list, cwd + '/probes_merged/' + ortho_file + '_merged.fasta', "fasta")

def merge_fasta_dir(in_dir, out_dir):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    for root, dirs, files in os.walk(in_dir):
        for exon_folder in dirs:
            seq_list = []
            exon_dir = os.path.join(root, exon_folder)
            #print exon_dir
            for exon_file in slistdir(exon_dir):
                exon_path = os.path.join(exon_dir, exon_file)
                #print exon_path
                for seq_record in SeqIO.parse(exon_path, "fasta"):
                    seq_list.append(seq_record)
            merged_fasta_name = out_dir + exon_folder + "_merged.fasta"
            SeqIO.write(seq_list, merged_fasta_name, "fasta")

#################################################################################################################

def multihold(commands):
    subprocess.call(commands, shell=True)

def parallel_exonerate():
    print "parallel_exonerate"
    commands = multirun_exonerate(refgenome)
    pool = mp.Pool(int(threads))
    pool.map(multihold, commands)

def parallel_find_in_alignment(exon_dir, nonref_dir, out_dir):
    print "parallel_find_in_alignment"
    commands = multirun_find_in_alignment(exon_dir, nonref_dir, out_dir)
    pool = mp.Pool(int(threads))
    pool.map(multihold, commands)

reftrantaxa = args.reference_taxon
refgenome = args.reference_genome
threads = args.threads
cwd = os.getcwd()

def split_mafft_alignments():
    print "split_mafft_alignments"
    locuscount = 1
    if os.path.exists('./locus_log.txt'):
        os.remove('./locus_log.txt')
    if not os.path.exists('split_mafft_fasta'):
        os.makedirs('split_mafft_fasta')
    for file in slistdir(orthofolder):
        for seq_record in SeqIO.parse(orthofolder + file, "fasta"):
            foldername = './split_mafft_fasta/' + 'Ortho' + str(locuscount)
            if not os.path.exists(foldername):
                os.makedirs(foldername)
            fastaname = foldername + '/' + 'Ortho' + str(locuscount) + '_' + seq_record.id.replace("|", "_") + '.fasta'
            namelog = 'Ortho' + str(locuscount) + ' ' + seq_record.id + ' ' + str(len(seq_record.seq.ungap('-')))
            with open('locus_log.txt', "a") as f:
                f.write(namelog + '\n')
            SeqIO.write(seq_record.upper(), fastaname, "fasta")
        locuscount += 1

#new_extract_transcriptome_refs(reftrantaxa, './split_mafft_fasta/', False)
#new_extract_transcriptome_refs(reftrantaxa, split_nonref_folder, ref_orthos_folder  )
def new_extract_transcriptome_refs(reftrantaxa, in_dir, out_dir):
    """Searches reference sequences for sequences with a fasta header containing an exact match for the reference taxa provided by user and copys those sequences to a separate folder"""
    print "extract_transcriptome_refs"
    seperator = "@" #character separating taxon and gene id in fasta headers. Agalma uses "@".
    reftrantaxa_list = reftrantaxa.split(' ') #buffer input, convert to list
    top_ortho_dict = {} #dict for holding lists of orthos
    best_match_dict = {}
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open('locus_log.txt') as f:
        for line in f:
            lines = line.strip().split(' ')
            if lines[0] not in top_ortho_dict:
                top_ortho_dict[lines[0]] = []
            next_entry = lines[1], lines[2]
            top_ortho_dict[lines[0]].append(next_entry)
    for reftaxa in reftrantaxa_list:
        for key in top_ortho_dict: #parse dict by orthoID
            ortho_list = top_ortho_dict[key]
            for ortho in ortho_list:
                ortho_taxon = ortho[0].split(seperator)[0]
                if ortho_taxon == reftaxa:
                    if key not in best_match_dict:
                        best_match_dict[key] = ortho[0]
    for ortho_folder in slistdir(in_dir):
        if ortho_folder in best_match_dict:
            cur_ortho_dir = in_dir + ortho_folder + '/'
            for ortho_file in slistdir(cur_ortho_dir):
                reference_found = []
                for seq_record in SeqIO.parse(cur_ortho_dir + ortho_file, "fasta"):
                        best_match = best_match_dict[ortho_folder]
                        if seq_record.id == best_match:
                            print "\tBest match for " + ortho_folder + ": " + seq_record.id
                            seq_record.seq = seq_record.seq.ungap("-")
                            fasta_name = ortho_folder + "_" + seq_record.id + ".fasta"
                            seq_record.id = "%s_%s" % (ortho_folder, seq_record.id)
                            SeqIO.write(seq_record.upper(), './ref_orthos/' + fasta_name, "fasta")
        else:
            print "\tNo representative ortholog found for " + ortho_folder
    #sys.exit("Done")

def collect_by_taxon(in_dir):
    seq_list = []
    seq_dict = {}
    if os.path.exists('exons_by_taxon'):
        shutil.rmtree('exons_by_taxon')
    os.makedirs('exons_by_taxon')
    for root, dirs, files in os.walk(in_dir, topdown=False):
        files = [f for f in files if not f[0] == '.']
        for name in files:
            in_fasta = os.path.join(root, name)
            seq_iterator = SeqIO.parse(in_fasta, "fasta")
            for record in seq_iterator:
                #record.id z.B.: Ortho1_Lionepha_casta@274393
                taxon_name = re.sub(r'Ortho[0-9]+_Exon[0-9]+_', '', record.id)
                taxon_name = re.sub(r'@[0-9]+$', '', taxon_name)
                with open('exons_by_taxon/' + taxon_name + ".fasta", "a") as f:
                    SeqIO.write(record, f, "fasta")

def walk_exon_ref_folders(in_dir, cdhit_file, out_dir):
    if os.path.exists("All_exons.fasta"):
        os.remove("All_exons.fasta")
    for root, dirs, files in os.walk(in_dir, topdown=False):
        files = [f for f in files if not f[0] == '.']
        for name in files:
            file_name = os.path.join(root, name)
            for seq_record in SeqIO.parse(file_name, "fasta"):
                with open("All_exons.fasta", "a+") as f:
                    SeqIO.write(seq_record, f, "fasta")
    cd_hit_exon_refs()
    split_cd_hit_out(cdhit_file, out_dir)
    os.rename(in_dir, "ref_exons_with_duplicates")
    os.rename(out_dir, in_dir)

def cd_hit_exon_refs():
    #subprocess.check_call(["cd-hit-est -i All_exons.fasta -o CDout.fasta  -c 0.9 -d 0 -M 16000 -n 8"], shell=True)
    cmd = "cd-hit-est -i All_exons.fasta -o CDout.fasta  -c 0.9 -d 0 -M 16000 -n 8"
    subprocess.check_call(shlex.split(cmd))

def split_cd_hit_out(cdhit_file, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for seq_record in SeqIO.parse(cdhit_file, "fasta"):
        print str(seq_record.id)
        names_split = str(seq_record.id).split("_")
        refexon_name = names_split[0].replace("Ref","")
        ortho_name = names_split[1]
        new_name = ortho_name + "_" + refexon_name
        ortho_dir = os.path.join(out_dir, ortho_name)
        if not os.path.exists(ortho_dir):
            os.makedirs(ortho_dir)
        new_dir = os.path.join(out_dir, ortho_name, new_name)
        with open(new_dir + "_Reference.fasta", "w+") as f:
                    SeqIO.write(seq_record, f, "fasta")

def append_intron_possition(in_dir, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    ortho_exon_dict = {}
    list_of_used_orthos = []
    for root, dirs, files in os.walk(in_dir, topdown=False):
        files = [f for f in files if not f[0] == '.']
        for name in files:
            file_name = os.path.join(root, name)
            ortho_exon = name.split("_merged")[0]
            ortho_name = ortho_exon.split("_")[0]
            exon_name =  ortho_exon.split("_")[1]
            if ortho_name not in ortho_exon_dict:
                print ortho_name + " not in dictionary"
                add_list = []
                add_list.append(file_name)
                ortho_exon_dict[ortho_name] = add_list
            else:
                print ortho_name + " found!"
                add_list = ortho_exon_dict[ortho_name]
                add_list.append(file_name)
                ortho_exon_dict[ortho_name] = add_list
    print ortho_exon_dict
    for key in ortho_exon_dict:
        ortho_list_for_rename = ortho_exon_dict[key]
        count_exons = len(ortho_list_for_rename)
        if count_exons != 1:
            for ortho in ortho_list_for_rename:
                new_ortho_name = ortho.split(root)[1]
                ortho_exon = new_ortho_name.split("_merged")[0]
                ortho_name = ortho_exon.split("_")[0]
                exon_name =  ortho_exon.split("_")[1]
                exon_number = int(exon_name.split("Exon")[1])
                if exon_number == 1 or exon_number == count_exons:
                    plus_intron_pos = "TI_" + ortho_exon + "_merged.fasta"
                    shutil.copy(ortho, out_dir + plus_intron_pos)
                else:
                    plus_intron_pos =  "II_" + ortho_exon + "_merged.fasta"
                    shutil.copy(ortho, out_dir + plus_intron_pos)
        else:
            ortho = ortho_list_for_rename[0]
            new_ortho_name = ortho.split(root)[1]
            ortho_exon = new_ortho_name.split("_merged")[0]
            plus_intron_pos =  "TT_" + ortho_exon + "_merged.fasta"
            shutil.copy(ortho, out_dir + plus_intron_pos)


def filter_by_taxa_and_length(in_dir, out_dir, min_taxa, min_length):
    print "filter_by_taxa_and_length"
    """Filters fasta files by minimum number of taxa and average sequence length, and adds this information to the file name"""
    if os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        shutil.rmtree(out_dir)
    if not os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        os.mkdir(out_dir)
    for root, dirs, files in os.walk(in_dir): #
        for name in files:     #This starts a loop of every file inside the directory
            if ".fasta" in name and "Ortho" in name: #this lets you capture only files you care about in the name list
                filename = os.path.join(root, name)
                count = 0
                count_list = [] #this is an empty list object
                sequences = []     #makes an empty list to be filled later
                for record in SeqIO.parse(filename, "fasta"):
                    length_seq = len(record.seq)
                    count_list.append(length_seq) #This fills the empty list we created above with the list of lengths for a file
                    count = count + 1
                    sequences.append(record)    #fills the empty list made above with all the records made in the command on the previous line
                round_mean = "%.2f" % round(np.mean(count_list),2)
                mean = np.mean(count_list)
                print "%s\n\tTaxa: %s Mean length:%s" % (name, str(len(sequences)), mean)
                if float(round_mean) < min_length:
                    print "\tMean length too short (%s)" % (round_mean)
                    #print "MEAN: %s Min length: %s" % (round_mean, min_length)
                if len(sequences) >= min_taxa and float(round_mean) >= min_length:
                    num_tax = len(count_list)
                    filename_fixed = filename.replace("./","")
                    output_name = out_dir + name #outfolder gets us to the right directory and filename determines what to name to new file
                    SeqIO.write(sequences, output_name, "fasta")    #takes objects from the list 'sequences' and write the relevant bits into a fasta file
                    new_name = "Len_" + str(round_mean) + "_" + "NumTax_" + str(num_tax) +"_" + name
                    os.rename(output_name, out_dir + new_name)
                    print "\tOrtho %s passed." % (name)
                else:
                    print "\tOrtho %s skipped!!!" % (name)

#align_with_macse("./probes_merged_Full_15_May_2017/")
def align_with_macse(in_dir, path_macse, out_dir):
    """Aligns merged fasta files using MACSE"""
    print "align_with_macse"
    commands = []
    if os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        shutil.rmtree(out_dir)
    if not os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        os.mkdir(out_dir)
    for fasta_file in slistdir(in_dir):
        fasta_file_path = os.path.join(in_dir, fasta_file)
        output_file_path = os.path.join(out_dir, fasta_file)
        command = 'java -jar -Xmx16000m %s -prog alignSequences -seq %s -out_NT %s -gc_def 1 -fs 10 -stop 15' % (path_macse, fasta_file_path, output_file_path)
        commands.append(command)
    return commands

def parallel_macse(macse_jar_path, in_dir, out_dir):
    print "parallel_macse"
    #cmd1 = "java -jar -Xmx16000m /Applications/macse_v1.2.jar -prog alignSequences -seq ./probes_merged_Full_15_May_2017/Len_996.80_NumTax_5_Ortho8262_Exon5_merged.fasta -out_NT ./macse_alignments/Len_996.80_NumTax_5_Ortho8262_Exon5_merged.fasta -gc_def 1 -fs 10 -stop 15"
    commands = align_with_macse(in_dir, macse_jar_path, out_dir)
    pool = mp.Pool(int(threads))
    pool.map(multihold, commands)

def clean_up_macse_alignments(in_dir):
    """Removes !s from MACSE alignments"""
    print "clean_up_macse_files"
    for file_name in slistdir(in_dir):
            if "merged.fasta" in file_name:
                in_file = os.path.join(in_dir, file_name)
                with open(in_file, "r+") as f:
                    for line in f:
                        new_line = line.replace("!","-")
                        new_file_name = in_file.replace("merged.fasta", "aligned.fasta")
                        with open(new_file_name, "a+") as g:
                            g.write(new_line)
                os.remove(in_file)

def align_with_mafft(in_dir, path_macse, out_dir):
    """Aligns merged fasta files using MAFFT"""
    print "align_with_mafft"
    commands = []
    if os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        shutil.rmtree(out_dir)
    if not os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        os.mkdir(out_dir)
    for fasta_file in slistdir(in_dir):
        fasta_file_path = os.path.join(in_dir, fasta_file)
        output_file_path = os.path.join(out_dir, fasta_file)
        #z.B.: "/usr/local/bin/mafft --nuc --thread 20 --ep 0 --genafpair --maxiterate 10000 /Volumes/Basement_Shelves/EverythingisTerrible/agalma/multalign-184/clusters/homologs_184_4012.fa > /Volumes/Basement_Shelves/EverythingisTerrible/agalma/analyses/Bembs_1/184/alignments/homologs_184_4012.fa"
        command = "%s --nuc --thread 1 --ep 0 --genafpair --maxiterate 10000 %s > %s" % (path_mafft, fasta_file_path, output_file_path)
        commands.append(command)
    return commands
    #print commands

def parallel_mafft(path_mafft, in_dir, out_dir):
    print "parallel_mafft"
    #cmd1 = "java -jar -Xmx16000m /Applications/macse_v1.2.jar -prog alignSequences -seq ./probes_merged_Full_15_May_2017/Len_996.80_NumTax_5_Ortho8262_Exon5_merged.fasta -out_NT ./macse_alignments/Len_996.80_NumTax_5_Ortho8262_Exon5_merged.fasta -gc_def 1 -fs 10 -stop 15"
    commands = align_with_mafft(in_dir, path_mafft, out_dir)
    pool = mp.Pool(int(threads))
    pool.map(multihold, commands)

# def clean_up_mafft_alignments(in_dir):
#     """Removes !s from MACSE alignments"""
#     print "clean_up_macse_files"
#     for file_name in slistdir(in_dir):
#             if "merged.fasta" in file_name:
#                 in_file = os.path.join(in_dir, file_name)
#                 with open(in_file, "r+") as f:
#                     for line in f:
#                         new_line = line.replace("!","-")
#                         new_file_name = in_file.replace("merged.fasta", "aligned.fasta")
#                         with open(new_file_name, "a+") as g:
#                             g.write(new_line)
#                 os.remove(in_file)


def get_longest_ORF(in_dir, out_dir, reference):
    print "get_longest_ORF"
    commands_lo = []
    commands_pd = []
    if os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        shutil.rmtree(out_dir)
    if not os.path.exists(out_dir): ##If the directory defined by out_folder doesn't exist, then make it (see next line), otherwise ignore it and print 'Exists!' (see else)
        os.mkdir(out_dir)
    if reference == True:
        for fasta_file in slistdir(in_dir):
            #fasta_path = os.path.join(in_dir, fasta_file)
            fasta_path = os.path.abspath( os.path.join(in_dir, fasta_file) )
            td_lo_command = "%s -t %s" % (transdecoder_LO_path, fasta_path)
            commands_lo.append(td_lo_command)
            td_pd_command = "%s -t %s  --no_refine_starts --single_best_only"    % (transdecoder_Pd_path, fasta_path)
            commands_pd.append(td_pd_command)

    else:
        for root, dirs, files in os.walk(in_dir, topdown=False):
            files = [f for f in files if not f[0] == '.']
            for name in files:
                # print root
                # print name
                fasta_path = os.path.abspath(os.path.join(root, name))
                td_lo_command = "%s -t %s" % (transdecoder_LO_path, fasta_path)
                commands_lo.append(td_lo_command)
                td_pd_command = "%s -t %s  --no_refine_starts --single_best_only"    % (transdecoder_Pd_path, fasta_path)
                commands_pd.append(td_pd_command)


    print commands_lo[0]
    print "\n"
    print commands_pd[0]
    print "\n****************************************************\n"
    return commands_lo, commands_pd

def parallel_transdecoder(in_dir, out_dir, reference=True):
    print "parallel_transdecoder"
    commands_lo, commands_pd = get_longest_ORF(in_dir, out_dir, reference )
    #sys.exit("Done")
    with cd(out_dir):
        pool1 = mp.Pool(int(threads))
        pool1.map(multihold, commands_lo)
        pool2 = mp.Pool(int(threads))
        pool2.map(multihold, commands_pd)

#clean_up_orfs(orf_folder , ref_orthos_folder, orf_log_name)
def clean_up_orfs(orf_dir, ortho_dir, orf_log_name, reference=True):
    print "clean_up_orfs"
    old_ortho_dir = ortho_dir.replace("./","").replace("/","") + "_transcripts"
    if os.path.exists(orf_log_name):
        os.remove(orf_log_name)
    print old_ortho_dir
    if os.path.exists(old_ortho_dir):
         shutil.rmtree(old_ortho_dir)
    os.rename(ortho_dir, old_ortho_dir)
    if not os.path.exists(ortho_dir):
            os.mkdir(ortho_dir)
    if reference == True:
        for orf_file in slistdir(orf_dir):
            if ".cds" in orf_file:
                orf_path = os.path.join(orf_dir, orf_file)
                for seq_record in SeqIO.parse(orf_path, "fasta"):
                    header_list = seq_record.id.split("::") #Parse TransDecoder header
                    seq_record.id = header_list[0].replace(">", "")
                    ortho_path = os.path.join(ortho_dir, seq_record.id)
                    print ortho_path
                    description = seq_record.description
                    seq_record.description = "" #Remove unnecessary description
                    SeqIO.write(seq_record, ortho_path, "fasta")
                    orf_type = description.split("type:")[1]
                    print orf_type
                    orf_type = "%s\tORF_Type: %s" % (seq_record.id, orf_type.split(" len:")[0])
                    with open(orf_log_name, "a") as f:
                        f.write(orf_type + '\n')
    else:
        for orf_file in slistdir(orf_dir):
            if ".cds" in orf_file:
                orf_path = os.path.join(orf_dir, orf_file)

                ortho_id = orf_file.split("_")[0]
                ortho_dir_path = os.path.join(ortho_dir, ortho_id)
                ortho_path = os.path.join(ortho_dir_path, orf_file)
                if not os.path.exists( ortho_dir_path ):
                    os.mkdir(ortho_dir_path)
                    print "Making ortho dir for %s" % (ortho_dir_path)
                for seq_record in SeqIO.parse(orf_path, "fasta"):

					# Old format:
					# >Gene.1::Ortho3_Lionepha_casta.fasta::g.1::m.1 Gene.1::Ortho3_Lionepha_casta.fasta::g.1  ORF type:complete len:413 (-) Ortho3_Lionepha_casta.fasta:647-1885(-)
					# New format:
					# >Ortho2_Bembidion_testatum.fasta.p1 GENE.Ortho2_Bembidion_testatum.fasta~~Ortho2_Bembidion_testatum.fasta.p1  ORF type:complete len:267 (+),score=336.49 Ortho2_Bembidion_testatum.fasta:187-987(+)


                    header_list = seq_record.id.split(" ") #Parse TransDecoder header
                    seq_record.id = header_list[0].replace(">", "")
                    seq_record.description = ""
                    ortho_full_path = os.path.join(ortho_dir, ortho_id, seq_record.id  )
                    print "\tWriting file: " + ortho_full_path
                    SeqIO.write(seq_record, ortho_full_path, "fasta")

    # Rename old ortho directory
    # Make new ortho directory (with the same name as the old directory)
    # Parse file header
    #    >Gene.1::Ortho6_Bembidion_sp_nr_transversale@47770::g.1::m.1 Gene.1::Ortho6_Bembidion_sp_nr_transversale@47770::g.1
    #    ORF type:complete len:377 (-) Ortho6_Bembidion_sp_nr_transversale@47770:275-1405(-)
    # Write new fasta file
    # [Log ORF type and length to text file]


#Nqc_walk_exon_ref_folders("./ref_exons/" , "All_orfs.fasta" , "CDout_orf.fasta", "./ref_exons_undup/")
def run_cd_hit(in_dir, combined_fasta_name , cdhit_file, out_dir):
    if os.path.exists(combined_fasta_name):
        os.remove(combined_fasta_name)
    for root, dirs, files in os.walk(in_dir, topdown=False):
        files = [f for f in files if not f[0] == '.']
        for name in files:
            file_name = os.path.join(root, name)
            for seq_record in SeqIO.parse(file_name, "fasta"):
                seq_record.seq = seq_record.seq.ungap("-")
                with open(combined_fasta_name, "a+") as f:
                    SeqIO.write(seq_record, f, "fasta")
    cd_hit_call( combined_fasta_name, cdhit_file )
    split_cd_hit_out(cdhit_file, out_dir)
    dup_dir_name = in_dir.replace("./","").replace("/","") + "_with_duplicates"

    if os.path.exists(dup_dir_name):
         shutil.rmtree(dup_dir_name)

    os.rename(in_dir, dup_dir_name)

    os.rename(out_dir, in_dir)

def split_cd_hit_out(cdhit_file, out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for seq_record in SeqIO.parse(cdhit_file, "fasta"):
        print str(seq_record.id)
        ortho_dir = os.path.join(out_dir, seq_record.id)
        with open(ortho_dir + ".fasta", "w+") as f:
                    SeqIO.write(seq_record, f, "fasta")

def cd_hit_call(in_fasta, out_fasta):
    #subprocess.check_call(["cd-hit-est -i All_exons.fasta -o CDout.fasta  -c 0.9 -d 0 -M 16000 -n 8"], shell=True)
    cmd = "cd-hit-est -i %s -o %s  -c 0.9 -d 0 -M 16000 -n 8" % (in_fasta, out_fasta)
    subprocess.check_call(shlex.split(cmd))




def main():
    start = time.time()
    main_ortho_folder = orthofolder

    #FASTA folders
    ref_orthos_folder = "./ref_orthos/"
    split_nonref_folder = "./split_mafft_fasta/"
    merged_probe_folder = "./probes_merged/"
    filtered_probe_folder = "./probes_filtered/"
    ref_exon_folder = "./ref_exons/"
    ref_exon_folder_position = "./ref_exons_position"
    probe_seq_folder = "./probe_seqs/"
    out_dir_macse = "./macse_alignments/"
    orf_folder = "./orfs/"
    #Temp folders
    undupped_exons = "./ref_exons_undup/"
    undupped_orthos = "./ref_orthos_undup/"
    #Exonerate outfile folders
    ref_exonerate_folder = "./t2g_exn_files/"
    paralogy_check_dir = "./exon2genome_exn_files"
    probe_exonerate_folder = "./probe_exn_files/"
    #Log files
    orf_log_name = "TransDecoder_log.txt"
    #Fasta files
    merged_exons_name = "All_exons.fasta"
    ref_exon_file = "./all_ref_exons.fasta"
    cdhit_exon_name = "CDout_exon.fasta"
    merged_orfs_name = "All_orfs.fasta"
    cdhit_orf_name = "CDout_orf.fasta"
    nonref_orf_folder = "./nonref_orfs/"
    nonref_log_name = "TransDecoder_NonRef_log.txt"
    #nonref_ortho_folder = "./nonref_orthos/"

    #clean_up_working(ref_exonerate_folder)

    #Preparing alignments
    split_mafft_alignments()
    new_extract_transcriptome_refs(reftrantaxa, split_nonref_folder , ref_orthos_folder  )
    ungap_split(split_nonref_folder)


    parallel_transdecoder(split_nonref_folder, nonref_orf_folder, False )
    clean_up_orfs(nonref_orf_folder , split_nonref_folder, nonref_log_name, False)

    #sys.exit("Done")

    #Converting transcripts to ORFs
    parallel_transdecoder(ref_orthos_folder, orf_folder )
    clean_up_orfs(orf_folder , ref_orthos_folder, orf_log_name)

    #Remove duplicate ORFs
    run_cd_hit(ref_orthos_folder, merged_orfs_name , cdhit_orf_name, undupped_orthos)

    #Find exons
    parallel_server_exonerate(ref_orthos_folder, refgenome, ref_exonerate_folder, 10)
    split_exonerate_parse(ref_exonerate_folder, ref_exon_folder, 'Exon', 0.95)
    ungap_split(split_nonref_folder)

    #Remove paralogs
    server_run_exonerate_paralogs(ref_exon_folder, refgenome, paralogy_check_dir, 6)

    remove_paralogs(paralogy_check_dir, ref_exon_folder, 'Exon', 0.76, 40)

    #sys.exit("Done")

    #Remove duplicate exons
    run_cd_hit(ref_exon_folder, merged_exons_name , cdhit_exon_name, undupped_exons)

    #Find exons in other sequences
    parallel_find_in_alignment(ref_exon_folder, split_nonref_folder, probe_exonerate_folder)
    qc_exon_to_fasta( probe_exonerate_folder, probe_seq_folder, 0.6 )

    #Organize output
    merge_fasta_dir(probe_seq_folder, merged_probe_folder)
    collect_by_taxon(merged_probe_folder)
    filter_by_taxa_and_length(merged_probe_folder, filtered_probe_folder, 3, 120)

    #Align exons
    #parallel_macse(path_macse, merged_probe_folder, out_dir_macse) #UNTESTED
    #macse_alignments_folder = "./macse_alignments_test_folder/"
    #clean_up_macse_alignments(macse_alignments_folder)


if __name__ == "__main__":
    main()
