# Exon_machine

Generates a list of candidate probe loci from a set of aligned or unaligned transcriptome sequences in a specified folder using a reference genome to determine exon bondaries. Transcriptome sequences should be provided as individual multi-fasta files, one file per ortholog, in a folder in the current working directory.

### Dependencies: 
  Exonerate, Transdecoder, CD-HIT-EST, Biopython

### Running test data

From the program folder, run:

`./ExonMachine.py -T 4 -t Bembidion_sp_nr_transversale -g test_genome.fasta -o test_orthologs/`

### Arguments

  -reference_taxon REFERENCE_TAXON, -t REFERENCE_TAXON
                        Specify which taxon in your ortholog clusters to use
                        as the source for reference sequences.
 
 -reference_genome REFERENCE_GENOME, -g REFERENCE_GENOME
                        Specify a fasta file to use as the reference genome.
  
  -ortho_directory ORTHO_DIRECTORY, -o ORTHO_DIRECTORY
                        Directory containing ortholog fasta files.

  -T THREADS, --threads THREADS
                        Number of threads to use for exonerate multithreading
                        (default=1).
