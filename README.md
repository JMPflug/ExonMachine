# Exon_machine

Generates a list of candidate probe loci from a set of aligned or unaligned transcriptome sequences in a specified folder using a reference genome to determine exon bondaries. Transcriptome sequences should be provided as individual multi-fasta files, one file per ortholog, in a folder in the current working directory.

### Dependencies: 

Exonerate, Transdecoder, CD-HIT-EST, Biopython

Exonerate and CD-HIT-EST should be in your PATH. Biopython can be installed using pip `pip install biopython`. The absolute path to both Transdecoder scripts (LongOrfs and Predict) should be specified in the `Configuration` at the beginning of the ExonMachine script.

Exonerate: https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate

Transdecoder: https://github.com/TransDecoder/TransDecoder/wiki

CD-HIT-EST: https://github.com/weizhongli/cdhit

Biopython: https://biopython.org/wiki/Download

### Data Requirements: 

All sequences should be in FASTA format.

#### Genome

The genome should be in a single FASTA file.

#### Ortholog Clusters

All genes in a ortholog group should be in a single FASTA file with a unique name, and all ortholog groups should be in the same directory. Note: Each cluster will be assigned a new unique ID number by default.

FASTA headers should only contain the a unique name for each specimen. This name should be the same across all clusters. Each cluster must contain the reference specimen specified by `-t`. Other missing specimen are fine.


### Running test data

From the program folder, run:

`./ExonMachine.py -T 4 -t Bembidion_sp_nr_transversale -g test_genome.fasta -o test_orthologs/`


### Output

Apologies that the output of ExonMachine is ~~somewhat~~ extremely excessive (this will be fixed in later versions).

Individual FASTA files containing the probes for each specimen can be found in the `probe_seqs` folder. Probes grouped by locus can be found in `probes_merged`. Descriptions of the other files and folders will be added here soon; however, the majority of these are intermediate files that are mainly useful for debugging.

### Arguments

-reference_taxon REFERENCE_TAXON, -t REFERENCE_TAXON:
Specify which taxon in your ortholog clusters to use as the source for reference sequences.
 
-reference_genome REFERENCE_GENOME, -g REFERENCE_GENOME:
Specify a fasta file to use as the reference genome.
  
-ortho_directory ORTHO_DIRECTORY, -o ORTHO_DIRECTORY:
Directory containing ortholog fasta files.

-T THREADS, --threads THREADS:
Number of threads to use for exonerate multithreading (default=1).
