#!/usr/bin/env python

import sys
import re
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
import glob
import argparse

#TODO:
# Add helper function to deal with DS_store issue
# Clean up directory traverses with os.walk() and os.path.join()
#	Make a traverse function and replace all the hard-coded listdir calls

parser = argparse.ArgumentParser(description="Generates a list of candidate probe loci from a set of aligned or unaligned\
											transcriptome sequences in a specified folder using a reference genome to\
											determine exon bondaries. Transcriptome sequences should be provided as\
											individual multi-fasta files, one file per ortholog, in a folder in the current\
											working directory.\
											Dependencies: Exonerate, Biopython")\

 											
parser.add_argument("-reference_taxon", "-t", required=True, help="Specify which taxon in your ortholog clusters to use as the source for reference sequences.",
					type=str)

parser.add_argument("-reference_genome", "-g", required=True, help="Specify a fasta file to use as the reference genome.",
					type=str)

parser.add_argument("-out_directory", "-o", required=True, help="Name to give output directory",
					type=str)

parser.add_argument("-T", "--threads", default=1, help="Number of threads to use for exonerate multithreading (default=1).",
					type=int)

args = parser.parse_args()

reftrantaxa = args.reference_taxon
refgenome = args.reference_genome
threads = args.threads

cwd = os.getcwd()

def slistdir(directory):
    """A specialized version of os.listdir() that ignores files that
    start with a leading period."""
    filelist = os.listdir(directory)
    return [x for x in filelist
            if not (x.startswith('.'))]

def extract_reference_transcriptome_seqs(reftrantaxa, dir):
	"""Reads a directory of multi-fasta orthologs as Bio.Seq objects
	and writes them 
	"""
	print "extract_reference_transcriptome_seqs"
	orthocount = 1
	if not os.path.exists('ref_orthos'):
		os.makedirs('ref_orthos')
	for file in slistdir(dir):
		for seq_record in SeqIO.parse(dir + file, "fasta"):
			if reftrantaxa in seq_record.id:
				seq_record_unpiped = seq_record.id.replace("|", "_")
				name = 'Ortho' +  str(orthocount) + '_' + seq_record_unpiped + '.fasta'
				SeqIO.write(seq_record, './ref_orthos/' + name, "fasta")
				orthocount += 1

def run_exonerate(refgenome): #exonerate_genome_to_transcriptome
	"""Compares reference genome and reference transcriptome to locate exons in orthologs.
	"""
	print "run_exonerate"
	if not os.path.exists('./t2g_exn_files/'):
		os.makedirs('./t2g_exn_files/')
	for file in slistdir('./ref_orthos/'):
		if 'DS_Store' not in file:
			ffile = './ref_orthos/' + file
			command = ['exonerate', ffile, refgenome, '-n', '1', '-m', 'est2genome']
			outfile = './t2g_exn_files/' + file + '.txt'
			with open(outfile, "w+") as f:
				subprocess.call(command, stdout=f)

def first_exonerate_parse(dir, newdir, prefix):
	print "first_exonerate_parse"
	cwd = os.getcwd()
	if not os.path.exists(cwd + newdir):
		os.makedirs(cwd + newdir)
	if not os.path.exists(cwd + '/merged_exons/'):
		os.makedirs(cwd + '/merged_exons/')
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


def join_exons(dir, dest, name):
	print "join_exons"
	os.chdir(dir)
	newcwd = os.getcwd()
	if os.path.exists(name):
		os.remove(name)
	for ffile in os.listdir(newcwd):
		if "DS_Store" not in ffile:
			with open(dest + name, "a") as f:	
				with open(ffile, "r") as g:
						f.write(g.read())
	os.chdir(cwd)

 
def split_mafft_alignments():
	print "split_mafft_alignments"
	locuscount = 1
	if os.path.exists('./locus_log.txt'):
		os.remove('./locus_log.txt')
	if not os.path.exists('split_mafft_fasta'):
		os.makedirs('split_mafft_fasta')
	for file in os.listdir('./orthologs/'):
		if "DS_Store" not in file:
			for seq_record in SeqIO.parse('./orthologs/' + file, "fasta"):
				foldername = './split_mafft_fasta/' + 'Locus' + str(locuscount)
				if not os.path.exists(foldername):
					os.makedirs(foldername)
				fastaname = foldername + '/' + 'Ortho' + str(locuscount) + '_' + seq_record.id + '.fasta'
				namelog = 'Ortho' + str(locuscount) + ' ' + seq_record.id
				with open('locus_log.txt', "a") as f:
					f.write(namelog + '\n')
				SeqIO.write(seq_record, fastaname, "fasta")
			locuscount += 1

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
								newname = OFloc + '/' + locuslist[0] + '_' + ortho
								oldname = OFloc + '/' + ortho 
								print oldname
								print newname
								os.rename(oldname, newname)

def process_locuslog(reftrantaxa):
	print "process_locuslog"
	if os.path.exists('./Locus_log_RefTaxa.txt'):
		os.remove('./Locus_log_RefTaxa.txt')
	with open('locus_log.txt') as f:
		for line in f:
			if reftrantaxa in line:
				wl = line.strip()
				with open('Locus_log_RefTaxa.txt', "a") as g:
					g.write(wl + '\n')

def ungap_split():
	print "ungap_split"
	for file in os.listdir('./split_mafft_fasta'):
		if "Locus" in file:
			for locusfolder in os.listdir('./split_mafft_fasta/' + file):
				for seq_record in SeqIO.parse('./split_mafft_fasta/' + file + '/' + locusfolder, "fasta"):
					newseq = seq_record.seq.ungap("-")
					os.remove('./split_mafft_fasta/' + file + '/' + locusfolder)
					with open('./split_mafft_fasta/' + file + '/' + locusfolder, "w+") as f:
						f.write('>' + locusfolder + '\n')
						f.write(str(newseq))

def find_in_alignment():
	print "find_in_alignment"
	if not os.path.exists('probes'):
		os.makedirs('probes')
	for file in os.listdir('./split_mafft_fasta/'):
		if 'Locus' in file:
			for subfile in os.listdir('./split_mafft_fasta/' + file):
				for ffile in os.listdir('./ref_exons/'):
					if "DS_Store" not in ffile:
						for fffile in os.listdir('./ref_exons/'  + ffile):
							if 'Locus' in fffile:
								fn = fffile.split('_')
								sn = subfile.split('_')
								if fn[0] == sn[0]:
									drffile =  './ref_exons/'  + ffile + '/' + fffile
									drsubfile = './split_mafft_fasta/' + file + '/' + subfile
									ryo = ">%ti %td\\n%tas\\n"
									command = ['exonerate', drffile, drsubfile, '-n', '1', '-m', 'affine:bestfit', '-E', '--ryo', ryo]
									ffs = fffile.split('_')
									newdir = './probes/' + file + '_' + ffs[1]
									if not os.path.exists(newdir):
										os.makedirs(newdir)
									outfile = newdir + '/' + fffile + '_vs_' + subfile 
									if os.path.exists(outfile):
										os.remove(outfile)
									with open(outfile, "w+") as f:
										subprocess.call(command, stdout=f)

def exon_to_fasta():
	print "exon_to_fasta"
	cwd = os.getcwd()
	if not os.path.exists('./probe_seqs'):
		os.makedirs('./probe_seqs')
	for file in os.listdir('./probes'):
		if 'Exon' in file:
			print file
			for subfile in os.listdir('./probes/' + file):
				#print file + " " + subfile
				if 'Ortho' in subfile:
					EL = './probes/' + file + '/'
					os.chdir(cwd + '/probes/' + file)
					result = SearchIO.parse(subfile, 'exonerate-text')
					for h in result:
						hit =  h[0][0][0]
						hitseq =  hit.hit
						fastaname = 'Probe_' + subfile
						if not os.path.exists(cwd + '/probe_seqs/' + file):
							os.makedirs(cwd + '/probe_seqs/' + file)
						SeqIO.write(hitseq, cwd + '/probe_seqs/' + file + '/' + fastaname, "fasta")
					os.chdir(cwd)

def ungap_split2(rootdir):
	print "ungap_split2"
	for file in os.listdir(rootdir):
		if "DS_Store" not in file:
			for locusfolder in os.listdir(rootdir + file):
				for seq_record in SeqIO.parse(rootdir + file + '/' + locusfolder, "fasta"):
					newseq = seq_record.seq.ungap("-")
					os.remove(rootdir + file + '/' + locusfolder)
					with open(rootdir + file + '/' + locusfolder, "w+") as f:
						f.write('>' + locusfolder + '\n')
						f.write(str(newseq))

def multirun_exonerate(refgenome):
	print "multirun_exonerate"
	if not os.path.exists('./t2g_exn_files/'):
		os.makedirs('./t2g_exn_files/')
	commands = []
	for file in os.listdir('./ref_orthos/'):
		if 'DS_Store' not in file:
			ffile = './ref_orthos/' + file
			outfile = './t2g_exn_files/' + file + '.txt'
			command = 'exonerate ' + "'" + ffile + "' " + refgenome + ' -n ' + '1 ' + '-m ' + 'est2genome ' + '> ' + "'" + outfile + "'"
			commands.append(command)
			# with open(outfile, "w+") as f:
			# 	subprocess.call(command, stdout=f)
	return commands


def multihold(commands):
	print "multihold"
	print "WARNING: Using Ctrl-c to stop this program will not work due to a bug in python multiprocessing. Use Ctrl-z and kill instead."
	subprocess.call(commands, shell=True)

def parallel_exonerate():
	print "parallel_exonerate"
	commands = multirun_exonerate(refgenome)
	pool = mp.Pool(int(threads))
	pool.map(multihold, commands)

#Check here for problem
def multirun_find_in_alignment():
	print "multirun_find_in_alignment"
	if not os.path.exists('probes'):
		os.makedirs('probes')
	commands = []
	for file in os.listdir('./split_mafft_fasta/'):
		#print file
		if 'Locus' in file:
			for subfile in os.listdir('./split_mafft_fasta/' + file):
				for ffile in os.listdir('./ref_exons/'):
					if "DS_Store" not in ffile:
						for fffile in os.listdir('./ref_exons/'  + ffile):
							if 'RefExon' in fffile:
								#print str(fffile) + "			fffile"
								#print str(subfile) + "			subfile"
								fn = fffile.split('_')
								sn = subfile.split('_')
								if fn[1] == sn[0]:
									drffile =  './ref_exons/'  + ffile + '/' + fffile
									drsubfile = './split_mafft_fasta/' + file + '/' + subfile
									ffs = fffile.split('_')
									print ffs
									print fffile
									newdir = './probes/' + file + '_' + ffs[0].replace("Ref","")
									if not os.path.exists(newdir):
										os.makedirs(newdir)
									outfile = newdir + '/' + fffile + '_vs_' + subfile 
									if os.path.exists(outfile):
										os.remove(outfile)
									command = 'exonerate ' + " '" + drffile + "' '" + drsubfile + "'" + ' -n 1 -m affine:bestfit -E > ' + "'" + outfile + "'"
									commands.append(command)
									# with open(outfile, "w+") as f:
									# 	subprocess.call(command, stdout=f)
	return commands

def parallel_find_in_alignment():
	print "parallel_find_in_alignment"
	commands = multirun_find_in_alignment()
	#print str(commands) + "Commands"
	pool = mp.Pool(int(threads))
	pool.map(multihold, commands)


def main():

	extract_reference_transcriptome_seqs(reftrantaxa, './orthologs/')
	commands = multirun_exonerate(refgenome)
	#parallel_exonerate()
	first_exonerate_parse('/t2g_exn_files/', '/ref_exons/', 'RefExon')

	split_mafft_alignments()
	rename_exons()
	process_locuslog(reftrantaxa)
	ungap_split()
	parallel_find_in_alignment()
	exon_to_fasta()
	ungap_split2('./probe_seqs/')

if __name__ == "__main__":
	main()



